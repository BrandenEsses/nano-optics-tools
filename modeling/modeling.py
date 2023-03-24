from flask import request, session, Blueprint, render_template
import matlab.engine
import numpy as np
import plotly.graph_objs as go
import plotly.offline as offline
from config import MODELING_MATLAB_FILES_PATH

modeling_blueprint = Blueprint('modeling', __name__, url_prefix='/modeling', template_folder="templates")

@modeling_blueprint.route('/')
def modeling():
    try:
        minT = session["minT"]
        maxT = session["maxT"]
        stepsT = session["stepsT"]
        minwn = session["minwn"]
        maxwn = session["maxwn"]
        stepswn = session["stepswn"]
        demod = session["demod"]
    except:
        minT = 300
        maxT = 1000
        stepsT = 7
        minwn = 800
        maxwn = 1600
        stepswn = 400
        demod = 2
    return render_template("modeling.html", figures=False, title="Modeling",
                           minT=minT, maxT=maxT, stepsT=stepsT,minwn=minwn,maxwn=maxwn,stepswn=stepswn, demod=demod)


@modeling_blueprint.route('/plot', methods=['POST'])
def plot():
    eng = matlab.engine.start_matlab()
    eng.cd(MODELING_MATLAB_FILES_PATH)
    minT = int(request.form.get("minT"))
    maxT = int(request.form.get("maxT"))
    stepsT = int(request.form.get("stepsT"))
    minwn = int(request.form.get("minwn"))
    maxwn = int(request.form.get("maxwn"))
    stepswn = int(request.form.get("stepswn"))
    demod = int(request.form.get("demod"))
    session["minT"] = minT
    session["maxT"] = maxT
    session["stepsT"] = stepsT
    session["minwn"] = minwn
    session["maxwn"] = maxwn
    session["stepswn"] = stepswn
    session["demod"] = demod
    T = np.linspace(minT, maxT, stepsT)
    nu =np.linspace(minwn, maxwn, stepswn)
    matlab_T = matlab.double(T)
    matlab_nu = matlab.double(nu)
    matlab_demod = matlab.double(demod)
    [snbar, amplitude, phase] = eng.Spheroid4Raschke(matlab_nu, matlab_T, matlab_demod, nargout=3)
    phase = np.array(phase)
    snbar = np.array(snbar)

    #Plot snbar at each T
    snbar_amp_data = []
    for i in range(len(matlab_T)):
        snbar_amp_data.append(go.Scatter(x=nu, y=np.abs(snbar[:, i]), mode='lines', name=str(round(T[i],2)) + "K"))
    layout = go.Layout(title=r'$\huge\textrm{Simualted s-SNOM Signal Amplitude}$',
                       xaxis={'title': r'$\Large\textrm{Wavenumber (cm-1)}$', 'automargin': True},
                       yaxis={'title': r'$\Large\textrm{Intensity (AU)}$', 'automargin': True},height=750, width=1000,
                       template="none")
    figure = go.Figure(data=snbar_amp_data, layout=layout)
    figure.update_layout(
        font=dict(
            family="Arial",
            size=22,  # Set the font size here
            color="Black"
        )
    )
    config = {'responsive': True}
    plot1_div = offline.plot(figure, auto_open=False, output_type='div', config=config)
    # Plot the phase at each T
    phase_data = []
    for i in range(len(T)):
        f = 1 - (T[i] - 300) * 3e-5;
        phase_data.append(go.Scatter(x=nu, y=phase[:,i], mode='lines', name=str(round(T[i],2))+"K; f=" + str(round(f,4))))
    phase_data.append(go.Scatter(x=nu, y=phase[:, 0]-phase[:, -1], mode='lines', name=str(round(T[0],2)) + "K" + " - " +str(round(T[-1],2)) + "K"))
    layout = go.Layout(title=r'$\huge\textrm{Simualted s-SNOM Signal Phase}$', xaxis={'title': r'$\Large\textrm{Wavenumber (cm-1)}$', 'automargin':True}, yaxis={'title': r'$\Large\textrm{Angle (Rad)}$', 'automargin':True},template="none",height=750,width=1000)
    figure = go.Figure(data=phase_data, layout=layout)
    figure.update_layout(
        font=dict(
            family="Arial",
            size=22,  # Set the font size here
            color="Black"
        )
    )
    config = {'responsive': True}
    plot2_div = offline.plot(figure, auto_open=False, output_type='div', config=config)
    return render_template('modeling.html', plot1_div=plot1_div,plot2_div=plot2_div,figures=True, title="Modeling",
                           minT=minT, maxT=maxT, stepsT=stepsT,minwn=minwn,maxwn=maxwn,stepswn=stepswn, demod=demod)