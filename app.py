from flask import Flask, render_template
from config import SECRET_KEY, UPLOAD_FOLDER

app = Flask(__name__)
app.config['SECRET_KEY'] = SECRET_KEY
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
@app.route('/')
def index():  # put application's code here
    return render_template("index.html", title="Nano-Optics Tools")

from nanoftir.nanoftir import nanoftir_blueprint
app.register_blueprint(nanoftir_blueprint)
from calculators.calculators import calculators_blueprint
app.register_blueprint(calculators_blueprint)
from pumpprobe.pumpprobe import pumpprobe_blueprint
app.register_blueprint(pumpprobe_blueprint)
from modeling.modeling import modeling_blueprint
app.register_blueprint(modeling_blueprint)
from conversions.conversions import conversions_blueprint
app.register_blueprint(conversions_blueprint)
from ape.ape import ape_blueprint
app.register_blueprint(ape_blueprint)

if __name__ == '__main__':
    app.run()
