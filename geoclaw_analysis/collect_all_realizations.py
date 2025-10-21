"""
Make html file of all realizations
slip, dtopo, hmax together

TODO:  dtopo_plot.png and subfaults_plot.png must be moved from
each geoclaw_output run directory to the png_files directory.

"""
import os

png_files = os.path.abspath('../png_files')


html_file = 'all_realizations_with_hmax.html'
html = open(html_file, 'w')
html.write('<html>\n')
html.write('<h1>All realizations</h1>')
html.write("""
<ul>
<li> <a href="#Mw8.6">Mw 8.6 realizations</a>
<li> <a href="#Mw8.8">Mw 8.8 realizations</a>
<li> <a href="#Mw9.0">Mw 9.0 realizations</a>
<li> <a href="#Mw9.2">Mw 9.2 realizations</a>
</ul>
<p>
""")

for Mw in [8.6,8.8,9.0,9.2]:
    html.write('<div id=Mw%s>' % Mw)
    html.write('<h2>Mw %s realizations</h2>' % Mw)
    for runno in range(100):
        html.write('<b>Mw %s, Realization %s</b><p>\n' % (Mw,runno))
        html.write('<img src="fine_%s/run_%s/subfaults_plot.png" height=400>&nbsp;\n' \
            % (Mw,runno))
        html.write('<img src="fine_%s/run_%s/dtopo_plot.png" height=400>\n' \
            % (Mw,runno))
        html.write('<img src="fine_%s/run_%s/run%s_Mw%s_fg3.png" height=400><p>\n' \
            % (Mw,runno,runno,Mw,))

html.write('</html>\n')
html.close()

print ("Created %s" % html_file)
print ("Move this file to the png_files directory to view")

