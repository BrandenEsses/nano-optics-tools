from flask import Flask, render_template, request, redirect, session, Blueprint
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import plotly.graph_objs as go
import plotly.offline as offline

ape_blueprint = Blueprint('ape_fitting', __name__, url_prefix='/ape_fitting', template_folder="templates")

def model_func(x, A, B, C, D, E):
    return E * x**4 + D * x**3 + C * x**2 + B * x + A

@ape_blueprint.route('/', methods=['GET','POST'])
def ape_fitting():
    return render_template("ape.html", title="APE Fitting", figures = False)

@ape_blueprint.route('/plot', methods=['POST'])
def ape_plotting():
    file = request.files['ape_file']
    file_content = np.loadtxt(file, delimiter = ",")
    wavelength = file_content[2:,4]
    powerdiode = file_content[2:,3]
    params, params_covariance = curve_fit(model_func, wavelength, powerdiode, p0=None)
    # Plotting the original data points
    # plt.scatter(wavelength, powerdiode, label='Data')

    # Plotting the fitted curve
    fit_start = np.min(wavelength)
    fit_end = np.max(wavelength)
    x_range = np.linspace(fit_start, fit_end, 100)
    y_pred = model_func(x_range, params[0], params[1], params[2], params[3], params[4])
    # plt.plot(x_range, y_pred, 'r-', label=)
    data = [go.Scatter(x=wavelength, y=powerdiode, mode='markers', name="Data"),
            go.Scatter(x=x_range, y=y_pred, mode='lines', marker=dict(color="Red"),name="Fit")]
    layout = go.Layout(title=r'$\large A + Bx + Cx^2 + Dx^3 + Ex^4$',
                       yaxis={'title': r'$\Large\textrm{Powerdiode}$', 'automargin': True},
                       xaxis={'title': r'$\Large\textrm{Wavelength (nm)}$', 'automargin': True}, height=750,
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
    figure.add_annotation(
        text=('A=%5.14f, B=%5.14f,\n C=%5.14f, D=%5.14f, E=%5.14f' % tuple(params))
        , showarrow=False
        , x=0.05
        , y=1
        , xref='paper'
        , yref='paper'
        , xanchor='left'
        , yanchor='top'
        , font=dict(size=12, color="black")
        , align="left")
    plot_div = offline.plot(figure, auto_open=False, output_type='div')
    return render_template("ape.html",title="APE Fitting", figures=True, plot_div=plot_div)

