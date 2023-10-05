from flask import Flask, render_template, request, redirect, session, Blueprint
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import plotly.graph_objs as go
import plotly.offline as offline

ape_blueprint = Blueprint('ape_fitting', __name__, url_prefix='/ape_fitting', template_folder="templates")

def func2(x, A, B):
    return B * x + A

def func3(x, A, B, C):
    return C * x**2 + B * x + A

def func4(x, A, B, C, D):
    return  D * x**3 + C * x**2 + B * x + A

def func5(x, A, B, C, D, E):
    return E * x**4 + D * x**3 + C * x**2 + B * x + A

@ape_blueprint.route('/', methods=['GET','POST'])
def ape_fitting():
    return render_template("ape.html", title="APE Fitting", figures = False)

@ape_blueprint.route('/plot', methods=['POST'])
def ape_plotting():
    plot_data = []
    column = request.form.get('data_select')
    file = request.files['ape_file']
    file_content = np.loadtxt(file, delimiter = ",", skiprows=1)
    file.seek(0)
    header = file.readline()
    string_header = header.decode(encoding='utf-8').replace('\r\n','').split(',')
    data_index = string_header.index(column)
    wavelength = file_content[:,0]
    selected_data = file_content[:,data_index]
    x_range = np.linspace(np.min(wavelength), np.max(wavelength), 100)
    plot_data.append(go.Scatter(x=wavelength, y=selected_data, mode='markers', name=column))

    popt2, pcov2 = curve_fit(func2, wavelength, selected_data, p0=None)
    popt3, pcov3 = curve_fit(func3, wavelength, selected_data, p0=None)
    popt4, pcov4 = curve_fit(func4, wavelength, selected_data, p0=None)
    popt5, pcov5 = curve_fit(func5, wavelength, selected_data, p0=None)

    # 2 parameters
    residuals2 = selected_data - func2(wavelength, *popt2)
    ss_res2 = np.sum(residuals2 ** 2)
    ss_tot2 = np.sum((selected_data - np.mean(selected_data)) ** 2)
    r_squared2 = 1 - (ss_res2 / ss_tot2)
    # print(r_squared2)
    y_pred2 = func2(x_range, popt2[0], popt2[1])
    plot_data.append(go.Scatter(x=x_range, y=y_pred2, mode='lines', name=r"$x$"))

    # 3 parameters
    residuals3 = selected_data - func3(wavelength, *popt3)
    ss_res3 = np.sum(residuals3 ** 2)
    ss_tot3 = np.sum((selected_data - np.mean(selected_data)) ** 2)
    r_squared3 = 1 - (ss_res3 / ss_tot3)
    y_pred3 = func3(x_range, popt3[0], popt3[1], popt3[2])
    plot_data.append(go.Scatter(x=x_range, y=y_pred3, mode='lines', name=r"$x^2$"))
    # print(r_squared3)

    # 4 parameters
    residuals4 = selected_data - func4(wavelength, *popt4)
    ss_res4 = np.sum(residuals4 ** 2)
    ss_tot4 = np.sum((selected_data - np.mean(selected_data)) ** 2)
    r_squared4 = 1 - (ss_res4 / ss_tot4)
    y_pred4 = func4(x_range, popt4[0], popt4[1], popt4[2], popt4[3])
    plot_data.append(go.Scatter(x=x_range, y=y_pred4, mode='lines', name=r"$x^3$"))
    # print(r_squared3)

    # 5 parameters
    residuals5 = selected_data - func5(wavelength, *popt5)
    ss_res5 = np.sum(residuals5 ** 2)
    ss_tot5 = np.sum((selected_data - np.mean(selected_data)) ** 2)
    r_squared5 = 1 - (ss_res5 / ss_tot5)
    y_pred5 = func5(x_range, popt5[0], popt5[1], popt5[2], popt5[3], popt5[4])
    plot_data.append(go.Scatter(x=x_range, y=y_pred5, mode='lines', name=r"$x^4$"))
    # print(r_squared3)

    layout = go.Layout(title=r'$\Large\textrm{Fitting}$',
                       yaxis={'title': r'$\textrm{' + column + '}$', 'automargin': True, 'tickformat':'none'},
                       xaxis={'title': r'$\textrm{Wavelength (cm}^{-1}\textrm{)}$', 'automargin': True, 'tickformat':'none'},
                       template="none")

    figure = go.Figure(data=plot_data, layout=layout)
    figure.update_layout(
        font=dict(
            family="Arial",
            size=22,  # Set the font size here
            color="Black"
        )
    )
    A2, B2 = popt2
    A3, B3, C3 = popt3
    A4, B4, C4, D4 = popt4
    A5, B5, C5, D5, E5 = popt5
    plot_div = offline.plot(figure, auto_open=False, output_type='div')
    table_layout = go.Layout(title=r'$\textrm{Fit function: }A + Bx + Cx^2 + Dx^3 + Ex^4$')
    table = go.Figure(data=[go.Table(header=dict(values=['Order',f'$A$', f'$B$', f'$C$', f'$D$', f'$E$', f'$R^2$']),
                                   cells=dict(values=[['Linear','Quadratic','Cubic','Quartic'],[A2,A3,A4,A5], [B2,B3,B4,B5], ['-',C3,C4,C5], ['-','-',D4,D5],['-','-','-',E5],[r_squared2,r_squared3,r_squared4,r_squared5]]))], layout=table_layout)
    table_div = offline.plot(table, auto_open=False, output_type='div')
    return render_template("ape.html",title="APE Fitting", figures=True, plot_div=plot_div, table_div=table_div)

