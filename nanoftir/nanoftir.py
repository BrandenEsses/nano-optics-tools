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
import plotly.graph_objs as go
import plotly.offline as offline
import pandas as pd

nanoftir_blueprint = Blueprint('nanoftir', __name__, url_prefix='/nanoftir', template_folder="templates")

@nanoftir_blueprint.route('/plot', methods=['POST'])
def plot():
    if 'file' not in request.files:
        return 'No file uploaded'
    file = request.files['file']
    df = pd.read_csv(file,sep='\t',names=["stage_pos_mm", "amplitude"])
    df['stage_pos_mm'] = df['stage_pos_mm'].astype(float)
    df['amplitude'] = df['amplitude'].astype(float)

    # First plot just the interferogram
    data = [go.Scatter(x=df['stage_pos_mm'], y=df['amplitude'], mode='markers')]
    layout = go.Layout(title=r'$\huge\textrm{Interferogram}$', xaxis={'title': r'$\Large\textrm{Stage Position (mm)}$', 'automargin':True}, yaxis={'title': r'$\Large\textrm{Amplitude (V)}$', 'automargin':True},height=750,width=1000,template="none")
    figure = go.Figure(data=data, layout=layout)
    figure.update_layout(
        font=dict(
            family="Arial",
            size=22,  # Set the font size here
            color="Black"
        )
    )
    plot1_div = offline.plot(figure, auto_open=False, output_type='div')

    # Take FFT of data
    xvals = np.array(df['stage_pos_mm'])/10**4
    yvals = np.array(df['amplitude'])
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
    x_min_index = np.argmin(np.abs(xf - 500)) # Find index where 500 wavenumbers occurs
    max_index = np.argmax(yf)
    data = [go.Scatter(x=xf[x_min_index:], y=yf[x_min_index:], mode='markers')]
    layout = go.Layout(title=r'$\huge\textrm{Spectrum}$', xaxis={'title': r'$\Large\textrm{Wavenumbers }(\textrm{cm}^{-1})$', 'automargin':True, 'tickformat':'digits'}, yaxis={'title': r'$\Large\textrm{Intensity (AU)}$','automargin':True},height=750,width=1000, template="none")
    figure = go.Figure(data=data, layout=layout)
    figure.update_layout(
        font=dict(
            family="Arial",
            size=22,  # Set the font size here
            color="Black"
        )
    )
    figure.add_annotation(x=xf[max_index], y=yf[max_index], text=r"$\textrm{" + str(round(xf[max_index])) + " cm}^{-1}$")
    plot2_div = offline.plot(figure, auto_open=False, output_type='div')
    return render_template('nanoftir.html', plot1_div=plot1_div, plot2_div=plot2_div,figures=True, title="NanoFTIR")

@nanoftir_blueprint.route('/',methods=['GET', 'POST'])
def nanoftir():
    return render_template("nanoftir.html",title="NanoFTIR",figures=False)