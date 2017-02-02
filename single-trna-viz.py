import requests
from flask import Flask, render_template, request, session
import sqlite3 as sql
import os
import subprocess
import re
app = Flask(__name__)

def increment_count():
  con = sql.connect("sessions.db")
  cur = con.cursor()
  cur.execute('INSERT INTO sessions DEFAULT VALUES')
  con.commit()
  con.close()
  return(get_session_count())

def get_session_count():
  con = sql.connect("sessions.db")
  cur = con.cursor()
  result = cur.execute('SELECT * FROM sessions ORDER BY ID DESC LIMIT 1')
  session_count = 100000 + result.fetchall()[0][0]
  return(session_count)

@app.route("/single-trna-heatmap", methods=['GET', 'POST'])
def display_page():
  errors = []
  single_plot_path = ''
  paired_plot_path = ''
  if request.method == "POST":
    session_id = increment_count()
    try:
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
      shell_cmd = 'single-trna-heatmap/single-trna-heatmap.R {} {} {} {}'.format(input_seq, input_species_clade, input_isotype, session_id)

      single_plot_path = 'single-plot-' + str(session_id) + '.png'
      paired_plot_path = 'paired-plot-' + str(session_id) + '.png'

      print('Parsed args, now running...')
      if subprocess.call(shell_cmd, shell=True) == 0:
        subprocess.call('mv {} static'.format(single_plot_path), shell=True)
        subprocess.call('mv {} static'.format(paired_plot_path), shell=True)
      else:
        errors.append("Something went wrong.")
    except:
      errors.append("Something went wrong.")
  
  return render_template('trna-viz.html', errors=errors, single_plot=single_plot_path, paired_plot=paired_plot_path)

if __name__ == "__main__":
  app.run(host='0.0.0.0', port=5900)
