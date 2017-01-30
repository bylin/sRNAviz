import requests
from flask import Flask, render_template, request
import subprocess
import re
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
      input_seq = input_seq.strip()
      if re.search(input_seq, "[^agctAGCT]"): 
        errors.append("Input sequence must contain only [AGCTagct]")
        return render_template('trna-viz.html', errors=errors, results=results)

      input_species_clade = ''
      if 'inputSpecies' in request.form: input_species_clade = request.form['inputSpecies']
      if 'inputClade' in request.form: input_species_clade = request.form['inputClade']
      if len(input_species_clade) == 0:
        errors.append("Did not specify species or clade")
        return render_template('trna-viz.html', errors=errors, results=results)

      input_isotype = request.form.getlist('inputIsotype') # Flask receives this as multiple values with the same key (inputIsotype => Ala, inputIsotype => Cys, etc.)
      input_isotype = ','.join(input_isotype) # so, we have to format and pass to R script
      shell_cmd = 'single-trna-heatmap/single-trna-heatmap.R {} {} {}'.format(input_seq, input_species_clade, input_isotype)
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
