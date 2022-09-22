import scipy.signal
from flask import Flask, render_template, request, redirect, session, Blueprint
from matplotlib import pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq
from scipy.special import erf
import math
import os, csv
import io
import random
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from app import app

nanoftir_blueprint = Blueprint('nanoftir', __name__, url_prefix='/nanoftir', template_folder="templates")

file = None

def get_data():
    global file
    stream = io.StringIO(file.decode("UTF8"), newline=None)
    tsv_file = csv.reader(stream, delimiter="\t")
    xvals = []
    yvals = []
    for line in tsv_file:
       xvals.append(10*float(line[0])/10**4) # Extra factor of 10 due to LabView VI issue
       yvals.append(float(line[1]))
    normalxvals = xvals
    normalyvals = yvals
    yvals = yvals - np.mean(yvals) # Normalize data
    length_scale = xvals[1] - xvals[0]
    num_data = len(yvals)
    data_length = np.log2(len(yvals))
    is_power_two = np.ceil(data_length) == np.floor(data_length)
    if not is_power_two:
        num_zeros_add = 2**(round(np.log2(num_data)) + 1) - num_data
        left_append = np.zeros(np.floor(num_zeros_add/2).astype(int))
        right_append = np.zeros(np.ceil(num_zeros_add/2).astype(int))
        yvals = np.insert(yvals,0,left_append)
        yvals = np.append(yvals,right_append)
    xvals = np.arange(len(yvals))*length_scale - len(left_append)*length_scale
    max_square_index = np.argmax(yvals**2)
    max_square_position = xvals[max_square_index]
    max_square_value = (yvals[max_square_index])**2
    half_max_square_position = xvals[np.argmin(abs(yvals**2 - 0.5*max_square_value))]
    stddev = np.abs(half_max_square_position - max_square_position)
    error_func = 0.5*(erf((xvals - half_max_square_position + stddev)/(np.sqrt(2)*stddev))+1)
    exp_tc = np.amax(xvals - np.amin(xvals))/4
    exponential = np.exp(-1*(xvals - max_square_position)/exp_tc)
    pure_apodization = error_func * exponential
    apodization = error_func*exponential*yvals
    N = len(xvals)
    yf = np.abs(fft(yvals)[0:N // 2])
    xf = fftfreq(N,length_scale)[0:N // 2]
    data = (xf, yf, normalxvals, normalyvals)
    return data

@nanoftir_blueprint.route('/intfg.png')
def plot_png():
    fig = create_figure()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

def create_figure():
    xf,yf,normalxvals,normalyvals = get_data()
    xf = xf/2
    fig = Figure()
    axis1 = fig.add_subplot(1, 2, 2)
    axis1.set_title("Spectrum")
    axis1.set_xlabel("Wavenumber (cm-1)")
    axis1.set_ylabel("Intensity (au)")
    ymax = max(yf)
    peak_indices,peak_heights_dict = scipy.signal.find_peaks(yf,height=0.00001*ymax, distance=50)
    peak_heights = peak_heights_dict['peak_heights']
    peak_indices = peak_indices[np.flip(np.argsort(peak_heights))][:3]
    colors = ['r^','y^','g^','b^','m^','c^','k^']
    i = 0
    for peak in peak_indices:
        text = str(round(xf[peak])) + " cm-1"
        try:
            axis1.plot(xf[peak], yf[peak],colors[i], label=text)
        except:
            continue
        i += 1
    axis1.plot(xf, yf)
    axis1.legend()
    axis2 = fig.add_subplot(1, 2, 1)
    axis2.set_title("Interferogram")
    axis2.set_xlabel("Position (cm)")
    axis2.set_ylabel("Intensity (V)")
    axis2.ticklabel_format(axis="y", style="sci", scilimits=(-6, -6))
    axis2.plot(normalxvals, normalyvals)
    fig.set_dpi(150)
    fig.set_size_inches(8.8, 4.95)
    fig.tight_layout(pad=2.0)
    return fig

@nanoftir_blueprint.route('/',methods=['GET', 'POST'])
def nanoftir():  # put application's code here
    if request.method == 'POST':
        global file
        file = request.files['file'].stream.read()
        return render_template("nanoftir.html", img_src=True, title="NanoFTIR")
    return render_template("nanoftir.html", img_src=False, title="NanoFTIR")