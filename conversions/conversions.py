from flask import request, session, Blueprint, render_template, Response
import gzip
import xml.etree.ElementTree as ET
import base64
import struct
import numpy as np

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
            stage_pos_mm = np.linspace(sweep_start, sweep_end, len(sweep_data_decoded))
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
    return render_template("conversions.html", title="File Conversions")

@conversions_blueprint.route('/axz', methods=['POST'])
def convert_axz():
    axz_file = request.files["axz_file"]
    spectra, heightmaps = get_spectra_interferograms(axz_file)
    return Response(str(spectra) + '\n' + str(heightmaps),
        mimetype='text/plain',
        headers={'Content-disposition': 'attachment; filename=spectra_heightmaps.txt'})