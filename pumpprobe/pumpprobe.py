from flask import Flask, render_template, request, redirect, session, Blueprint, send_file
import numpy as np
import os, csv
import io
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy.optimize import curve_fit
import math
import plotly.graph_objs as go
import plotly.offline as offline
import pandas as pd

pumpprobe_blueprint = Blueprint('pumpprobe', __name__, url_prefix='/pumpprobe', template_folder="templates")

def exponential(x, a, b, c):
    return a * np.exp(-x/b) + c

def biexponential(x, a, b, c, d, e):
    return a * np.exp(-x/b) + c * np.exp(-x/d) + e

def triexponential(x, a, b, c, d, e, f, g):
    return a * np.exp(-x/b) + c * np.exp(-x/d) + e * np.exp(-x/f) + g

@pumpprobe_blueprint.route('/',methods=['GET', 'POST'])
def pumpprobe():  # put application's code here
    return render_template("pumpprobe.html", figures=False, title="Pump-Probe")

@pumpprobe_blueprint.route('/plot', methods=['POST'])
def plot():
    remove_from_front = int(request.form.get("removefront"))
    remove_from_end = int(request.form.get("removeend"))
    num_terms = int(request.form.get("fitnum"))
    files = request.files.getlist("file")
    dfx = pd.read_csv(files[0],sep='\t',names=["stage_pos_mm", "amplitude","1","2","3","4","5"])
    dfy = pd.read_csv(files[1], sep='\t', names=["stage_pos_mm", "amplitude","1","2","3","4","5"])
    xvalues = np.array(dfx["amplitude"])/10**6
    yvalues = np.array(dfy["amplitude"])/10**6
    stage_positions_mm = np.array(dfx["stage_pos_mm"])
    length_scale = stage_positions_mm[1] - stage_positions_mm[0]
    ps_spacing = length_scale * 2 / 2.99e11 * 1e12
    remove_from_front_index = round(remove_from_front / ps_spacing)
    max_cut = math.floor(len(stage_positions_mm) / 2)
    remove_from_end_index = len(stage_positions_mm) - round(remove_from_end / ps_spacing)
    xvalues = np.flip(xvalues)[remove_from_front_index:remove_from_end_index]
    yvalues = np.flip(yvalues)[remove_from_front_index:remove_from_end_index]
    rvalues = np.sqrt(np.array(xvalues) ** 2 + np.array(yvalues) ** 2)
    normalized_r = np.array(abs(rvalues - 0.5 * max(rvalues)))
    max_index = np.argmax(rvalues)
    max_index = np.argmax(rvalues) + round(0.1 * abs(max_index - len(xvalues)))
    min_index = max_index + np.argmin(rvalues[max_index:])
    t0_index = np.argmin(normalized_r[:max_index])
    t0_position = stage_positions_mm[t0_index]
    stage_positions_ps = ((np.array(stage_positions_mm) - t0_position) * 2 / 2.99e11 * 1e12)[
                         remove_from_front_index:remove_from_end_index]
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

    if num_terms == 1:
        a = popt[0]
        b = popt[1]
        c = popt[2]
        fit_string = str(round(a, 2)) + "*exp(-t/" + str(round(b, 2)) + ") + " + str(round(c, 2))
        if (a,b,c) == (0,1,0):
            fit_string = "Fit failure"
        data = [go.Scatter(x=stage_positions_ps, y=rvalues, mode='markers', name="Experimental Data"),go.Scatter(x=fit_x, y=[exponential(fit_x, a, b, c)], mode='markers',name="Fit")]
        layout = go.Layout(title=r'$\huge\textrm{Fitted Amplitude (single exponential)}$',
                           xaxis={'title': r'$\Large\textrm{Stage Position (ps)}$', 'automargin': True},
                           yaxis={'title': r'$\Large\textrm{Amplitude (V)}$', 'automargin': True}, height=750,
                           width=1000,
                           template="none")
        figure = go.Figure(data=data, layout=layout)
        figure.update_layout(
            font=dict(
                family="Arial",
                size=22,  # Set the font size here
                color="Black"
            )
        )
        figure.add_annotation(x=np.argmax(rvalues), y=1.05*max(rvalues),
                              text=fit_string, showarrow=False)
        plot2_div = offline.plot(figure, auto_open=False, output_type='div')
    if num_terms == 2:
        a = popt[0]
        b = popt[1]
        c = popt[2]
        d = popt[3]
        e = popt[4]
        fit_string = str(round(a,2)) + "*exp(-t/" + str(round(b,2)) + ") + " + str(round(c,2)) + "*exp(-t/" + str(round(d,2)) + ") + " + str(round(e,2))
        if (a,b,c,d,e) == (0,1,0,1,0):
            fit_string = "Fit failure"
        data = [go.Scatter(x=stage_positions_ps, y=rvalues, mode='markers', name="Experimental Data"),go.Scatter(x=fit_x, y=[biexponential(fit_x,a,b,c,d,e)], mode='markers',name="Fit")]
        layout = go.Layout(title=r'$\huge\textrm{Fitted Amplitude (biexponential)}$',
                           xaxis={'title': r'$\Large\textrm{Stage Position (ps)}$', 'automargin': True},
                           yaxis={'title': r'$\Large\textrm{Amplitude (V)}$', 'automargin': True}, height=750,
                           width=1000,
                           template="none")
        figure = go.Figure(data=data, layout=layout)
        figure.update_layout(
            font=dict(
                family="Arial",
                size=22,  # Set the font size here
                color="Black"
            )
        )
        figure.add_annotation(x=np.argmax(rvalues), y=1.05*max(rvalues),
                              text=fit_string, showarrow=False)
        plot2_div = offline.plot(figure, auto_open=False, output_type='div')
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
        data = [go.Scatter(x=stage_positions_ps, y=rvalues, mode='markers', name="Experimental Data"),go.Scatter(x=fit_x, y=triexponential(fit_x,a,b,c,d,e,f,g), mode='markers', name="Fit")]
        layout = go.Layout(title=r'$\huge\textrm{Fitted Amplitude (triexponential)}$',
                           xaxis={'title': r'$\Large\textrm{Stage Position (ps)}$', 'automargin': True},
                           yaxis={'title': r'$\Large\textrm{Amplitude (V)}$', 'automargin': True}, height=750,
                           width=1000,
                           template="none")
        figure = go.Figure(data=data, layout=layout)
        figure.add_annotation(x=np.argmax(rvalues), y=1.05*max(rvalues),
                              text=fit_string, showarrow=False)
        figure.update_layout(
            font=dict(
                family="Arial",
                size=22,  # Set the font size here
                color="Black"
            )
        )
        plot2_div = offline.plot(figure, auto_open=False, output_type='div')

    # Plot the amplitude
    data = [go.Scatter(x=stage_positions_ps, y=rvalues, mode='markers')]
    layout = go.Layout(title=r'$\huge\textrm{Amplitude}$',
                       xaxis={'title': r'$\Large\textrm{Stage Position (ps)}$', 'automargin': True},
                       yaxis={'title': r'$\Large\textrm{Amplitude (V)}$', 'automargin': True}, height=750, width=1000,
                       template="none")
    figure = go.Figure(data=data, layout=layout)
    figure.update_layout(
        font=dict(
            family="Arial",
            size=22,  # Set the font size here
            color="Black"
        )
    )
    plot1_div = offline.plot(figure, auto_open=False, output_type='div')


    return render_template("pumpprobe.html", img_src=False, title="Pump-Probe", plot1_div = plot1_div,plot2_div = plot2_div,figures=True)