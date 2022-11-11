from flask import Flask, render_template, request, redirect, session, Blueprint
import numpy as np
import os, csv
import io
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy.optimize import curve_fit

pumpprobe_blueprint = Blueprint('pumpprobe', __name__, url_prefix='/pumpprobe', template_folder="templates")

x_file = None
y_file = None

def biexponential(x, a, b, c, d, e):
    return a * np.exp(-x/b) + c * np.exp(-x/d) + e

def get_data_pumpprobe():
    global x_file
    global y_file
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
    xvalues = np.flip(xvalues)
    yvalues = np.flip(yvalues)
    rvalues = np.sqrt(np.array(xvalues)**2 + np.array(yvalues)**2)
    normalized_r = np.array(abs(rvalues - 0.5*max(rvalues)))
    max_index = np.argmax(rvalues)
    t0_index = np.argmin(normalized_r[:max_index])
    t0_position = stage_positions_mm[t0_index]
    stage_positions_ps = (np.array(stage_positions_mm) - t0_position) * 2 / 2.99e11 * 1e12
    min_index = max_index + np.argmin(rvalues[max_index:])
    fit_y = np.array(rvalues[max_index:min_index])
    fit_x = np.array(stage_positions_ps[max_index:min_index])
    popt, pcov = curve_fit(biexponential, fit_x, fit_y, p0=(max(rvalues),10,max(rvalues),10,0))
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
    a = popt[0]
    b = popt[1]
    c = popt[2]
    d = popt[3]
    e = popt[4]
    fit_string = str(round(a,2)) + "*exp(-t/" + str(round(b,2)) + ") + " + str(round(c,2)) + "*exp(-t/" + str(round(d,2)) + ") + " + str(round(e,2))
    axis3.plot(fit_x, biexponential(fit_x,a,b,c,d,e), linewidth=0.8, label=fit_string)
    axis3.legend()
    fig.set_dpi(170)
    fig.set_size_inches(4.95,12)
    fig.tight_layout(pad=2.0)
    return fig

@pumpprobe_blueprint.route('/',methods=['GET', 'POST'])
def pumpprobe():  # put application's code here
    if request.method == 'POST':
        global x_file
        global y_file
        files = request.files.getlist("file")
        x_file = files[0].stream.read()
        y_file = files[1].stream.read()
        return render_template("pumpprobe.html", img_src=True, title="Pump-Probe")
    return render_template("pumpprobe.html", img_src=False, title="Pump-Probe")