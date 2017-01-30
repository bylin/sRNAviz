import requests
from flask import Flask, render_template, request
import subprocess
app = Flask(__name__)

@app.route("/single-trna-heatmap", methods=['GET', 'POST'])

# def submit():
#   error = None
#   if request.method == 'POST':
#     if valid

def display_page():
  errors = []
  results = []
  if request.method == "POST":
    try:
      subprocess.call('rm static/single-plot.png static/paired-plot.png', shell=True)
      input_seq = request.form['inputSeq']
      input_species = request.form['inputSpecies']
      input_isotype = request.form['inputIsotype']
      input_best_isotype = request.form['inputBestIsotype']
      shell_cmd = 'single-trna-heatmap/single-trna-heatmap.R {} {} {} {}'.format(input_seq, input_species, input_isotype, input_best_isotype)
      print('Parsed args, now running...')
      if subprocess.call(shell_cmd, shell=True) == 0:
        subprocess.call('mv single-plot.png static', shell=True)
        subprocess.call('mv paired-plot.png static', shell=True)
        results = True
      else:
        results = False
        errors.append("Something went wrong.")
    except:
      errors.append("Something went wrong.")
  
  return render_template('trna-viz.html', errors=errors, results=results)

if __name__ == "__main__":
  app.run(host='0.0.0.0', port=5900)
