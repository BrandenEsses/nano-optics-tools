from flask import Flask, render_template, request, redirect, session
from matplotlib import pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq
import math
import os, csv
import io
import random
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
app = Flask(__name__)
app.config['SECRET_KEY'] = "issacnewton"

file = None

def get_data():
    global file
    stream = io.StringIO(file.decode("UTF8"), newline=None)
    tsv_file = csv.reader(stream, delimiter="\t")
    xvals = []
    yvals = []
    for line in tsv_file:
       xvals.append(float(line[0]))
       yvals.append(float(line[1]))
    yvals = yvals - np.mean(yvals) # Normalize data
    length_scale = xvals[1] - xvals[0]
    normalxvals = xvals
    normalyvals = yvals
    num_data = len(yvals)
    data_length = np.log2(len(yvals))
    is_power_two = np.ceil(data_length) == np.floor(data_length)
    if not is_power_two:
        num_zeros_add = 2**(round(np.log2(num_data)) + 1) - num_data
        left_append = np.zeros(np.floor(num_zeros_add/2).astype(int))
        right_append = np.zeros(np.ceil(num_zeros_add/2).astype(int))
        yvals = np.insert(yvals,0,left_append)
        yvals = np.append(yvals,right_append)
    new_num_data = len(yvals)
    yvals_zeros = np.zeros(new_num_data-256)
    max_location = np.argmax(yvals)
    max_value = yvals[max_location]
    half_max_index = np.argmin(yvals - 0.5*max_value)
    min_index = max_location - 128
    max_index = max_location + 128
    slice = yvals[min_index:max_index]
    yvals_zeros = np.insert(yvals_zeros,min_index,slice)
    xvals = np.arange(len(yvals)) * length_scale
    stddev = np.abs(half_max_index - max_location) * length_scale
    gaussian = 1/np.sqrt(2*np.pi*stddev**2) * np.exp(-(xvals - max_location*length_scale)**2/(2*stddev**2))
    apodization = yvals_zeros*gaussian
    right_half= apodization[max_location:]
    left_half = np.flip(apodization[:max_location])
    rotated = np.append(right_half,left_half)
    yf = np.absolute(fft(rotated)[0:num_data//2])
    xf = fftfreq(num_data, length_scale)[:num_data//2]
    data = (xf, yf, normalxvals, normalyvals)
    return data

@app.route('/intfg.png')
def plot_png():
    fig = create_figure()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

def create_figure():
    xs,ys,xvals,normalyvals = get_data()
    fig = Figure()
    axis1 = fig.add_subplot(1, 2, 2)
    axis1.set_title("Spectrum")
    axis1.set_xlabel("Frequency (au)")
    axis1.set_ylabel("Intensity (au)")
    axis1.plot(xs, ys)
    axis2 = fig.add_subplot(1, 2, 1)
    axis2.set_title("Interferogram")
    axis2.set_xlabel("Position (au)")
    axis2.set_ylabel("Intensity (au)")
    axis2.plot(xvals, normalyvals)
    fig.set_dpi(150)
    #fig.set_size_inches(11.2,6.3)
    fig.set_size_inches(8.8, 4.95)
    fig.tight_layout(pad=2.0)
    return fig

@app.route('/')
def index():  # put application's code here
    return render_template("index.html", title="Nano-Optics Tools")

@app.route('/nanoftir',methods=['GET', 'POST'])
def nanoftir():  # put application's code here
    if request.method == 'POST':
        global file
        file = request.files['file'].stream.read()
        return render_template("nanoftir.html",img_src=True, title="NanoFTIR")
    return render_template("nanoftir.html",img_src=False,title="NanoFTIR")

if __name__ == '__main__':
    app.run()
