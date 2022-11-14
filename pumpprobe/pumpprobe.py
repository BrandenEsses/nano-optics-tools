from flask import Flask, render_template, request, redirect, session, Blueprint, send_file
import numpy as np
import os, csv
import io
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy.optimize import curve_fit
import math

pumpprobe_blueprint = Blueprint('pumpprobe', __name__, url_prefix='/pumpprobe', template_folder="templates")

x_file = None
y_file = None
num_terms = None
remove_from_front = None
remove_from_end = None
max_cut = None

def exponential(x, a, b, c):
    return a * np.exp(-x/b) + c

def biexponential(x, a, b, c, d, e):
    return a * np.exp(-x/b) + c * np.exp(-x/d) + e

def triexponential(x, a, b, c, d, e, f, g):
    return a * np.exp(-x/b) + c * np.exp(-x/d) + e * np.exp(-x/f) + g

def get_data_pumpprobe():
    global x_file
    global y_file
    global remove_from_front
    global remove_from_end
    global max_cut
    x_stream = io.StringIO(x_file.decode("UTF8"), newline=None)
    x_tsv_file = csv.reader(x_stream, delimiter="\t")
    y_stream = io.StringIO(y_file.decode("UTF8"), newline=None)
    y_tsv_file = csv.reader(y_stream, delimiter="\t")
    xvalues = []
    yvalues = []
    stage_positions_mm = []
    for line in x_tsv_file:
        stage_positions_mm.append(float(line[0]))
        xvalues.append(float(line[1]))
    for line in y_tsv_file:
        yvalues.append(float(line[1]))
    ps_spacing = (stage_positions_mm[1] - stage_positions_mm[0]) * 2 / 2.99e11 * 1e12
    remove_from_front_index = round(remove_from_front/ps_spacing)
    max_cut = math.floor(len(stage_positions_mm)/2)
    remove_from_end_index = len(stage_positions_mm) - round(remove_from_end/ps_spacing)
    xvalues = np.flip(xvalues)[remove_from_front_index:remove_from_end_index]
    yvalues = np.flip(yvalues)[remove_from_front_index:remove_from_end_index]
    rvalues = np.sqrt(np.array(xvalues)**2 + np.array(yvalues)**2)
    normalized_r = np.array(abs(rvalues - 0.5*max(rvalues)))
    max_index = np.argmax(rvalues)
    max_index = np.argmax(rvalues) + round(0.1*abs(max_index-len(xvalues)))
    min_index = max_index + np.argmin(rvalues[max_index:])
    t0_index = np.argmin(normalized_r[:max_index])
    t0_position = stage_positions_mm[t0_index]
    stage_positions_ps = ((np.array(stage_positions_mm) - t0_position) * 2 / 2.99e11 * 1e12)[remove_from_front_index:remove_from_end_index]
    fit_y = np.array(rvalues[max_index:min_index])
    fit_x = np.array(stage_positions_ps[max_index:min_index])
    try:
        if num_terms == 1:
            popt, pcov = curve_fit(exponential, fit_x, fit_y, p0=(max(rvalues), 10, 0))
        if num_terms == 2:
            popt, pcov = curve_fit(biexponential, fit_x, fit_y, p0=(max(rvalues),10,max(rvalues),10,0))
        if num_terms == 3:
            popt, pcov = curve_fit(triexponential, fit_x, fit_y, p0=(max(rvalues),10,max(rvalues),10,max(rvalues),10,0))
    except:
        popt = []
        for i in range(2*num_terms + 1):
            if i%2 == 0:
                popt.append(0)
            else:
                popt.append(1)
    return (stage_positions_ps, xvalues, yvalues, rvalues, fit_x, popt)

@pumpprobe_blueprint.route('/pp1.png')
def plot_png():
    fig = create_figure()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

def create_figure():
    stage_positions_ps,xvalues,yvalues,rvalues,fit_x,popt = get_data_pumpprobe()
    fig = Figure()
    axis1 = fig.add_subplot(4, 1, 1)
    axis1.set_title("X Values")
    axis1.set_xlabel("Pump Delay (ps)")
    axis1.set_ylabel("Intensity (µV)")
    axis1.plot(stage_positions_ps, xvalues,linewidth=0.8)
    axis2 = fig.add_subplot(4, 1, 2)
    axis2.set_title("Y Values")
    axis2.set_xlabel("Pump Delay (ps)")
    axis2.set_ylabel("Intensity (µV)")
    axis2.plot(stage_positions_ps, yvalues,linewidth=0.8)
    axis3 = fig.add_subplot(4, 1, 3)
    axis3.set_title("Amplitude")
    axis3.set_xlabel("Pump Delay (ps)")
    axis3.set_ylabel("Intensity (µV)")
    axis3.plot(stage_positions_ps, rvalues,linewidth=0.8)
    # axis3.set_xticks(np.arange(stage_positions_ps[0],stage_positions_ps[-1]+10,10))
    axis3 = fig.add_subplot(4, 1, 4)
    axis3.set_title("Fitted Amplitude")
    axis3.set_xlabel("Pump Delay (ps)")
    axis3.set_ylabel("Intensity (µV)")
    axis3.plot(stage_positions_ps, rvalues, linewidth=0.8, label="Data")
    if num_terms == 1:
        a = popt[0]
        b = popt[1]
        c = popt[2]
        fit_string = str(round(a, 2)) + "*exp(-t/" + str(round(b, 2)) + ") + " + str(round(c, 2))
        if (a,b,c) == (0,1,0):
            fit_string = "Fit failure"
        axis3.plot(fit_x, exponential(fit_x, a, b, c), linewidth=0.8, label=fit_string)
    if num_terms == 2:
        a = popt[0]
        b = popt[1]
        c = popt[2]
        d = popt[3]
        e = popt[4]
        fit_string = str(round(a,2)) + "*exp(-t/" + str(round(b,2)) + ") + " + str(round(c,2)) + "*exp(-t/" + str(round(d,2)) + ") + " + str(round(e,2))
        if (a,b,c,d,e) == (0,1,0,1,0):
            fit_string = "Fit failure"
        axis3.plot(fit_x, biexponential(fit_x,a,b,c,d,e), linewidth=0.8, label=fit_string)
    if num_terms == 3:
        a = popt[0]
        b = popt[1]
        c = popt[2]
        d = popt[3]
        e = popt[4]
        f = popt[5]
        g = popt[6]
        fit_string = str(round(a,2)) + "*exp(-t/" + str(round(b,2)) + ") + " + str(round(c,2)) + "*exp(-t/" + str(round(d,2)) + ") + " + str(round(e,2)) + "*exp(-t/" + str(round(f,2)) + ") + " + str(round(g,2))
        if (a,b,c,d,e,f,g) == (0,1,0,1,0,1,0):
            fit_string = "Fit failure"
        axis3.plot(fit_x, triexponential(fit_x,a,b,c,d,e,f,g), linewidth=0.8, label=fit_string)

    axis3.legend(loc='lower left')
    fig.set_dpi(150)
    fig.set_size_inches(7,20)
    fig.tight_layout(pad=2.0)
    return fig

@pumpprobe_blueprint.route('/',methods=['GET', 'POST'])
def pumpprobe():  # put application's code here
    global x_file
    global y_file
    global remove_from_front
    global remove_from_end
    global num_terms
    global max_cut

    remove_from_front = 0
    remove_from_end = 0
    num_terms = 2

    try:
        session["front"]
        session["back"]
        session["numterms"]
    except:
        session["front"] = remove_from_front
        session["back"] = remove_from_end
        session["num_terms"] = num_terms

    if request.method == 'POST':
        remove_from_front = int(request.form.get("removefront"))
        remove_from_end = int(request.form.get("removeend"))
        num_terms = int(request.form.get("fitnum"))
        files = request.files.getlist("file")
        x_file = files[0].stream.read()
        y_file = files[1].stream.read()
        session["front"] = remove_from_front
        session["back"] = remove_from_end
        session["num_terms"] = num_terms

        return render_template("pumpprobe.html", img_src=True, title="Pump-Probe", numterms=session["num_terms"], front=session["front"], back=session["back"])

    return render_template("pumpprobe.html", img_src=False, title="Pump-Probe",numterms=session["num_terms"], front=session["front"], back=session["back"])