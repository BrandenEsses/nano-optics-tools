from flask import request, session, Blueprint, render_template, Response
import gzip
import xml.etree.ElementTree as ET
import base64
import struct
import numpy as np
from .get_spectra import get_spectrum_from_interferogram
import time
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.offline as offline

conversions_blueprint = Blueprint('conversions', __name__, url_prefix='/conversions', template_folder="templates")

def axz_to_xml(filepath):
    f = gzip.open(filepath, 'rb')
    file_content = f.read()
    xml = file_content.decode('UTF-16')
    return xml


def xml_to_root(xml):
    root = ET.fromstring(xml)
    return root


def base64string_to_decimal(string):
    decoded_string = base64.b64decode(string)
    scan_data_grouped_raw = []
    for j in range(1, len(decoded_string)):
        if j % 4 == 0:
            binary_numbers = decoded_string[j - 4:j]
            scan_data_grouped_raw.append(binary_numbers)
    scan_data_grouped_decimal = []
    for num in scan_data_grouped_raw:
        scan_data_grouped_decimal.append(struct.unpack('f', num)[0])
    return scan_data_grouped_decimal


def decimal_string_to_scan(string, scan_dimension):
    final_scan_data = []
    for k in range(1, len(string)):
        if k % scan_dimension == 0:
            final_scan_data.append(string[k - scan_dimension:k])
    return np.array(final_scan_data)


def get_datatype(root, datatype):
    for element in root[:]:
        actual_datatype = element.tag.split('}')[1]
        if actual_datatype == datatype:
            return (element)
    return None


def get_interferograms(snom_spectra):
    all_spectra_list = []
    for spectra in snom_spectra:
        label = spectra[0].text
        interferograms = spectra[13]
        for channel_scan in interferograms:
            channel_scan_dict = {}
            sweep_start = float(channel_scan[0].text)
            sweep_end = float(channel_scan[1].text)
            sweep_channel = channel_scan[8].text
            sweep_data_encoded = channel_scan[9].text
            sweep_data_decoded = base64string_to_decimal(sweep_data_encoded)
            stage_pos_mm = np.array(np.linspace(sweep_start, sweep_end, len(sweep_data_decoded)))
            channel_scan_dict.update({'label': label})
            channel_scan_dict.update({'data_channel': sweep_channel})
            channel_scan_dict.update({'stage_pos_mm': stage_pos_mm})
            channel_scan_dict.update({'data': sweep_data_decoded})
            all_spectra_list.append(channel_scan_dict)
    return all_spectra_list


def get_heightmaps(heightmaps):
    all_heightmaps_list = []
    for heightmap in heightmaps:
        heightmap_dict = {}
        label = heightmap.attrib['Label']
        data_channel = heightmap.attrib['DataChannel']
        dim_x = int(heightmap[3][0].text)
        dim_y = int(heightmap[3][1].text)
        try:
            unit = heightmap[5].text + heightmap[4].text
        except:
            unit = heightmap[4].text
        encoded_data = heightmap[11].text
        decoded_data = base64string_to_decimal(encoded_data)
        data = decimal_string_to_scan(decoded_data, dim_x)
        heightmap_dict.update({'label': label})
        heightmap_dict.update({'data_channel': data_channel})
        heightmap_dict.update({'dimensions': (dim_x, dim_y)})
        heightmap_dict.update({'data': data})
        all_heightmaps_list.append(heightmap_dict)
    return all_heightmaps_list


def get_spectra_interferograms(filepath):
    xml = axz_to_xml(filepath)
    root = xml_to_root(xml)
    snom_spectra = get_datatype(root, "SNOMSpectra")
    heightmaps = get_datatype(root, "HeightMaps")
    all_heightmaps_list = get_heightmaps(heightmaps)
    all_spectra_list = get_interferograms(snom_spectra)
    return [all_spectra_list, all_heightmaps_list]

@conversions_blueprint.route('/', methods=['GET','POST'])
def conversions():
    try:
        sample_label = session["axz_sample_label"]
        data_channel = session["axz_data_channel"]
        minWN = int(session["axz_minWN"])
        maxWN = int(session["axz_maxWN"])
        cutoff = int(session["axz_cutoff"])
        resLB = int(session["axz_resLB"])
        resUB = int(session["axz_resUB"])
    except:
        sample_label = ""
        data_channel = ""
        minWN = 1650
        maxWN = 1800
        cutoff = 10
        resLB = 0
        resUB = 0
    params = (sample_label, data_channel, minWN, maxWN, cutoff, resLB, resUB)
    return render_template("conversions.html", title="File Conversions", figures = False, form_values = params)

@conversions_blueprint.route('/axz', methods=['POST'])
def convert_axz():
    axz_file = request.files["axz_file"]
    axz_ref = request.files["axz_ref"]
    sample_label = request.form.get("sample_label")
    data_channel = request.form.get("data_channel")
    minWN = int(request.form.get("minWN"))
    maxWN = int(request.form.get("maxWN"))
    cutoff = int(request.form.get("cutoff"))
    resLB = int(request.form.get("resLB"))
    resUB = int(request.form.get("resUB"))
    session["axz_sample_label"] = sample_label
    session["axz_data_channel"] = data_channel
    session["axz_minWN"] = minWN
    session["axz_maxWN"] = maxWN
    session["axz_cutoff"] = cutoff
    session["axz_resLB"] = resLB
    session["axz_resUB"] = resUB
    spectra, heightmaps = get_spectra_interferograms(axz_file)
    ref_data = np.loadtxt(axz_ref, delimiter='\t')
    params = (minWN,maxWN,cutoff,resLB,resUB)
    processed_spectra = []
    start = time.time()
    for spectrum in spectra:
        if spectrum['data_channel'] != data_channel and data_channel != "":
            continue
        if spectrum["label"] != sample_label and sample_label != "":
            continue
        complete_label = spectrum["label"] + " - " + spectrum['data_channel']
        print(complete_label)
        samp_intfgm = np.array(spectrum['data']).astype(float)
        stage_pos_mm = np.array(spectrum['stage_pos_mm']).astype(float)
        samp_data = np.transpose(np.stack((10**6*stage_pos_mm,samp_intfgm))) # Factor of 1000 to get mm -> µm
        samp_data_len_rows = len(samp_data)
        ref_data_len_rows = len(ref_data)
        samp_data_len_cols = len(samp_data[0,:])
        ref_data_len_cols = len(ref_data[0,:])
        diff_rows = samp_data_len_rows - ref_data_len_rows
        ref_cols = ref_data_len_cols
        # if samp_data_len_rows > ref_data_len_rows:
        #     ref_data = np.append(ref_data, np.zeros((diff_rows,ref_cols)) + np.mean(ref_data, axis=0), axis = 0)
        # print(len(ref_data),len(samp_data))
        sampArray, refAvg = get_spectrum_from_interferogram(samp_data, ref_data, params)
        processed_spectra.append([complete_label,sampArray, refAvg])
        progress = time.time()
        print("Elapsed time: " + str(progress - start) + " seconds")
    plot_divs = []
    spectrum = processed_spectra[0]
    fig1 = make_subplots(specs=[[{"secondary_y": True}]])
    fig1.add_trace(
        go.Scatter(
            x=np.abs(spectrum[2][:, 0]),
            y=np.angle(spectrum[2][:, 1]),
            mode='lines',
            name="Phase"
        ),
        secondary_y=True,
    )

    fig1.add_trace(
        go.Scatter(
            x=np.abs(spectrum[2][:, 0]),
            y=np.abs(spectrum[2][:, 1]),
            mode='lines',
            name="Amplitude"
        ),
        secondary_y=False,
    )

    fig1.update_layout(
        title_text=spectrum[0] + " - Reference Spectrum",
        height=750,
        width=1000,
        template="none",
        font=dict(
            family="Arial",
            size=22,  # Set the font size here
            color="Black"
        ),
    )

    fig1.update_xaxes(title_text='Wavenumbers (cm-1)')
    fig1.update_yaxes(title_text='Amplitude (V)', secondary_y=False)
    fig1.update_yaxes(title_text="Phase (Rad)", secondary_y=True)
    plot_div = offline.plot(fig1, auto_open=False, output_type='div')
    plot_divs.append(plot_div)

    for spectrum in processed_spectra:
        fig2 = make_subplots(specs=[[{"secondary_y": True}]])
        fig2.add_trace(
            go.Scatter(
                x=np.abs(spectrum[1][:, 0]),
                y=np.abs(spectrum[1][:, 1]),
                mode='lines',
                name="Amplitude"
            ),
            secondary_y=True,
        )
        fig2.add_trace(
            go.Scatter(
                x=np.abs(spectrum[1][:, 0]),
                y=np.angle(spectrum[1][:, 1]),
                mode='lines',
                name="Phase"
            ),
            secondary_y=False,
        )

        fig2.update_layout(
            title_text=spectrum[0] +  " - Referenced Sample Spectrum",
            height=750,
            width=1000,
            template="none",
            font=dict(
                family="Arial",
                size=22,  # Set the font size here
                color="Black"
            ),
        )

        fig2.update_xaxes(title_text='Wavenumbers (cm-1)')
        fig2.update_yaxes(title_text='Amplitude (V)', secondary_y=False)
        fig2.update_yaxes(title_text="Phase (Rad)", secondary_y=True)
        plot_div = offline.plot(fig2, auto_open=False, output_type='div')
        plot_divs.append(plot_div)
    sample_label = session["axz_sample_label"]
    data_channel = session["axz_data_channel"]
    minWN = session["axz_minWN"]
    maxWN = session["axz_maxWN"]
    cutoff = session["axz_cutoff"]
    resLB = session["axz_resLB"]
    resUB = session["axz_resUB"]
    params = (sample_label, data_channel, minWN, maxWN, cutoff, resLB, resUB)
    return render_template("conversions.html", title="File Conversions", figures=True, plots = plot_divs, form_values = params)