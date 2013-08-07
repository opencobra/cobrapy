from os.path import dirname, abspath, join, isfile, isdir
from warnings import warn
from urllib2 import urlopen
import json

maps_dir = join(abspath(dirname(__file__)), "")
map_cache_dir = join(maps_dir, "map_cache", "")
static_dir = join(maps_dir, "static", "")


d3_filepath = join(static_dir, 'd3.v3.js')
map_js_filepath = join(static_dir, 'visbio_map.js')
with open(join(static_dir, 'map.css')) as infile:
    style = infile.read().replace("\n", " ")
ipython_html = """
<style>
.overlay {
  fill: none;
  pointer-events: all;
}
</style>
<button onclick="download_map()">Download svg</button>
<div id="map"></div>"""



map_download_url = "http://zakandrewking.github.io/visbio/maps/"


class Map(object):
    """View metabolic map"""
    def __init__(self, map_file="e-coli-core", flux={}):
        if map_file.endswith(".json"):
            warn("Map file name should not include .json")
        # if the file is not present attempt to download
        map_filename = join(maps_dir, "map_cache", map_file + ".json")
        if not isfile(map_filename):
            map_not_cached = 'Map "%s" not in cache. Attempting download from %s' % \
                (map_file, map_download_url)
            warn(map_not_cached)
            from urllib2 import urlopen
            if not isdir(map_cache_dir):
                from os import mkdir
                mkdir(map_cache_dir)
            download = urlopen(map_download_url + map_file + ".json")
            # TODO catch HTTP Error and raise more descriptive exception
            with open(map_filename, "w") as outfile:
                outfile.write(download.read())
        with open(map_filename) as f:
            self.map_json = f.read()
        self.flux = flux

    def _repr_html_(self):
        from IPython.display import HTML
        with open(html_filepath, 'r') as f:
            html = f.read()
        javascript = self._assemble_javascript()
        javascript += """\nvisBioMap.visualizeit(d3.select("#map"), map_data, style, flux, null, null, null);"""
        javascript += """\nsvg = document.getElementsByTagName("svg")[0];"""
        return html + """<script type="text/Javascript">%s</script>""" % javascript


    def _assemble_javascript(self):
        with open(d3_filepath, 'r') as f:
            d3 = f.read()
        with open(map_js_filepath, 'r') as f:
            map_js = f.read()
        javascript = "\n".join([
            "var " + d3, map_js,
            "var map_data = %s;" % self.map_json,
            "var flux = %s;" % json.dumps(self.flux),
            'var style = "%s";' % style])
        return javascript

    def create_standalone_html(self, outfile=None):
        with open(join(static_dir, "standalone_map_template.html")) as infile:
            template = infile.read()
        template = template.replace("__ALL_JAVASCRIPT__", self._assemble_javascript())
        with open("/home/aebrahim/test2.html", "w") as outfile:
            outfile.write(template)
    
    
    def view_browser(self):
        import webbrowser
        # TODO finish


    def run_server(self):
        """start a tornado server to display the map"""
        None  # TODO implement
