var visBioMap = (function(d3) {
    var maps = {};
    maps.version = 0.1;

    maps.visualizeit = function(selection, data, style, flux, flux2, metabolites, metabolites2) {
        maps.defaults = {};
        maps.defaults.path_color = "rgb(80, 80, 80)";

        maps.style_variables = get_style_variables(style);

        var decimal_format = d3.format('.1f');
        var decimal_format_3 = d3.format('.3f');
        maps.has_flux = false;
        maps.has_flux_comparison = false;
        maps.has_metabolites = false;
        maps.has_metabolite_deviation = false;

        // parse the data objects and attach values to map objects
        if (flux) {
            maps.has_flux = true;
            data = parse_flux_1(data, flux);
            if (flux2) {
                maps.has_flux_comparison = true;
                data = parse_flux_2(data, flux2);
            }
        }
        if (metabolites) {
            maps.has_metabolites = true;
            data = parse_metabolites_1(data, metabolites);
            if (metabolites2) {
                maps.has_metabolite_deviation = true;
                data = parse_metabolites_2(data, metabolites2);
            }
        }

        var map_w = data.max_map_w, map_h = data.max_map_h;
        var width = parseFloat(selection.style("width")) - 20;
        var height = parseFloat(selection.style("height")) - 20;
        var factor = Math.min(width/map_w,height/map_h);
        console.log('map_w '+decimal_format(map_w));
        console.log('map_h '+decimal_format(map_h));
        console.log('scale factor '+decimal_format_3(factor));

        maps.scale = {};
        maps.scale.x = d3.scale.linear()
            .domain([0, map_w])
            .range([(width - map_w*factor)/2, map_w*factor + (width - map_w*factor)/2]);
        maps.scale.y = d3.scale.linear()
            .domain([0, map_h])
            .range([(height - map_h*factor)/2, map_h*factor + (height - map_h*factor)/2]);
        maps.scale.x_size = d3.scale.linear()
            .domain([0, map_w])
            .range([0, map_w*factor]);
        maps.scale.y_size = d3.scale.linear()
            .domain([0, map_h])
            .range([0, map_h*factor]);
        maps.scale.size = d3.scale.linear()
            .domain([0, 1])
            .range([0, factor]);
        maps.scale.flux = d3.scale.linear()
            .domain([0, 40])
            .range([6, 6]);
        maps.scale.flux_fill = d3.scale.linear()
            .domain([0, 40, 200])
            .range([1, 1, 1]);
        maps.scale.flux_color = d3.scale.linear()
            .domain([0, 1, 20, 50])
            .range(["rgb(200,200,200)", "rgb(150,150,255)", "blue", "red"]);
        maps.scale.metabolite_concentration = d3.scale.linear()
            .domain([0, 10])
            .range([15, 200]);
        maps.scale.metabolite_color = d3.scale.linear()
            .domain([0, 1.2])
            .range(["#FEF0D9", "#B30000"]);

        var svg = selection
            .attr("style", "width:"+width+"px;height:"+height+"px;margin:0px auto")// ;border:3px solid black;")
            .append("svg")
        // TODO: add correct svg attributes (see '/Users/zaking/Dropbox/lab/optSwap/paper-2-GAPD/old figs/fig5-theoretical-production/')
            .attr("width", width)
            .attr("height", height)
            .attr("xmlns", "http://www.w3.org/2000/svg")
            .append("g")
            .call(d3.behavior.zoom().scaleExtent([1, 15]).on("zoom", zoom))
            .append("g");
        maps.svg = svg;

        svg.append("style")
            .attr("type", "text/css")
            .text(style);

        var g = svg.append('g')
            .attr('id', 'markers');

        g.append("marker")
            .attr("id", "end-triangle-path-color")
            .attr("markerHeight", 2.1)
            .attr("markerUnits", "strokeWidth")
            .attr("markerWidth", 2.1)
            .attr("orient", "auto")
            .attr("refX", 0)
            .attr("refY", 6)
            .attr("viewBox", "0 0 12 12")
            .append("path")
            .attr("d", "M 0 0 L 12 6 L 0 12 z")
            .attr("fill", maps.defaults.path_color)
        // stroke inheritance not supported in SVG 1.*
            .attr("stroke", maps.defaults.path_color);

        g.append("marker")
            .attr("id", "start-triangle-path-color")
            .attr("markerHeight", 2.0)
            .attr("markerUnits", "strokeWidth")
            .attr("markerWidth", 2.0)
            .attr("orient", "auto")
            .attr("refX", 12)
            .attr("refY", 6)
            .attr("viewBox", "0 0 12 12")
            .append("path")
            .attr("d", "M 12 0 L 0 6 L 12 12 z")
            .attr("fill", maps.defaults.path_color)
            .attr("stroke", maps.defaults.path_color);

        svg.append("rect")
            .attr("class", "overlay")
            .attr("width", width)
            .attr("height", height)
            .attr("style", "stroke:black;fill:none;");

        svg.append("g")
            .attr("id", "membranes")
            .selectAll("rect")
            .data(data.membrane_rectangles)
            .enter().append("rect")
            .attr("class", function(d){ return d.class })
            .attr("width", function(d){ return maps.scale.x_size(d.width) })
            .attr("height", function(d){ return maps.scale.y_size(d.height) })
            .attr("transform", function(d){return "translate("+maps.scale.x(d.x)+","+maps.scale.y(d.y)+")";})
            .style("stroke-width", function(d) { return maps.scale.size(10) })
            .attr('rx', function(d){ return maps.scale.x_size(20); })
            .attr('ry', function(d){ return maps.scale.x_size(20); });

        if (data.hasOwnProperty("metabolite_circles")) {
            console.log('metabolite circles');
            svg.append("g")
                .attr("id", "metabolite-circles")
                .selectAll("circle")
                .data(data.metabolite_circles)
                .enter().append("circle")
                .attr("r", function (d) {
                    sc = maps.scale.metabolite_concentration; // TODO: make a better scale
                    if (d.metabolite_concentration) {
                        var s;
                        if (d.should_size) s = maps.scale.size(sc(d.metabolite_concentration));
                        else s = maps.scale.size(0);
                        return s;
                    } else if (maps.has_metabolites) {
                        return maps.scale.size(10);
                    } else {
                        return maps.scale.size(d.r);
                    }
                })
                .attr("style", function (d) {
                    sc = maps.scale.metabolite_color;
                    if (d.metabolite_concentration) {
                        var a;
                        if (d.should_color) a = "fill:"+sc(d.metabolite_concentration) + ";" +
                            "stroke:black;stroke-width:0.5;";
                        else a = "fill:none;stroke:black;stroke-width:0.5;";
                        return a;
                    }
                    else if (maps.has_metabolites) {
                        return "fill:grey;stroke:none;stroke-width:0.5;";
                    }
                    else { return ""; }
                })
                .attr("transform", function(d){return "translate("+maps.scale.x(d.cx)+","+maps.scale.y(d.cy)+")";});
            if (maps.has_metabolite_deviation) {
                append_deviation_arcs(data.metabolite_circles);
            }
        }
        else if (data.hasOwnProperty("metabolite_paths")) {
            if (maps.has_metabolites) { alert('metabolites do not render w simpheny maps'); }
            console.log('metabolite paths');
            svg.append("g")
                .attr("id", "metabolite-paths")
                .selectAll("path")
                .data(data.metabolite_paths)
                .enter().append("path")
                .attr("d", function(d) { return scale_path(d.d); })
                .style("fill", "rgb(224, 134, 91)")
                .style("stroke", "rgb(162, 69, 16)")
                .style("stroke-width", String(maps.scale.size(2))+"px");
        }

        append_reaction_paths(data.reaction_paths);

        svg.append("g")
            .attr("id", "reaction-labels")
            .selectAll("text")
            .data(data.reaction_labels)
            .enter().append("text")
            .text(function(d) {
                var t = d.text;
                if (maps.has_flux_comparison)
                    t += " ("+decimal_format(d.flux1)+"/"+decimal_format(d.flux2)+": "+decimal_format(d.flux)+")";
                else if (d.flux) t += " ("+decimal_format(d.flux)+")";
                else if (maps.has_flux) t += " (0)";
                return t;
            })
            .attr("text-anchor", "start")
            .attr("font-size", function(d) {
                var s;
                if (maps.style_variables.hasOwnProperty('reaction_label_size')) {
                    s = maps.style_variables['reaction_label_size'];
                }
                else { s = 15; }
                return maps.scale.size(s);
            })
        // .attr("style", function(d){ if(!d.flux) return "visibility:hidden;"; else return ""; })
            .attr("transform", function(d){return "translate("+maps.scale.x(d.x)+","+maps.scale.y(d.y)+")";});

        svg.append("g")
            .attr("id", "misc-labels")
            .selectAll("text")
            .data(data.misc_labels)
            .enter().append("text")
            .text(function(d) { return d.text; })
            .attr("font-size", maps.scale.size(60))
            .attr("transform", function(d){return "translate("+maps.scale.x(d.x)+","+maps.scale.y(d.y)+")";});

        svg.append("g")
            .attr("id", "metabolite-labels")
            .selectAll("text")
            .data(data.metabolite_labels)
            .enter().append("text")
            .text(function(d) {
                var t = d.text;
                if (isNaN(d.metabolite_concentration)) {}
                else if (has_metabolite_deviation) {
                    var a = (isNaN(d.metabolite_concentration) ? "-" : decimal_format(d.metabolite_concentration));
                    var b = (isNaN(d.metabolite_deviation) ? "-" : decimal_format(d.metabolite_deviation));
                    t += " ("+a+" \xB1 "+b+"%)";
                }
                else if (d.metabolite_concentration) {
                    var a = (isNaN(d.metabolite_concentration) ? "-" : decimal_format(d.metabolite_concentration));
                    t += " ("+a+")";
                }
                else if (has_metabolites) t += " (0)";
                return t;
            })
            .attr("font-size", function(d) {
                if (d.metabolite_concentration) return maps.scale.size(30);
                else if (maps.has_metabolites) return maps.scale.size(20);
                else return maps.scale.size(20);
            })
            .style("visibility","visible")
            .attr("transform", function(d){return "translate("+maps.scale.x(d.x)+","+maps.scale.y(d.y)+")";});


        function zoom() {
            svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
        }

    }

    function get_style_variables(style) {
	return [];

        // var r = new RegExp("/\*(a-zA-Z)+\*/");
        // var r = /\/\*([a-zA-Z_]+)=([0-9.]+)\*\//g;
        var r = /\/\*\s([a-zA-Z_]+)\s=\s([0-9.]+)\s\*\//g;
        matches = r.exec(style);        // only returns the first one...need findAll
        params = {};
        params[matches[1]] = matches[2];
        return params
    }

    function parse_flux_1(data, flux) {
        data.reaction_paths = data.reaction_paths.map( function(o) {
            // console.log(d3.keys(flux));
            if (o.id in flux) {
                o.flux = parseFloat(flux[o.id]);
            }
            // else { console.log(o.id) }
            return o;
        });
        data.reaction_labels = data.reaction_labels.map( function(o) {
            if (o.text in flux) {
                // TODO: make sure text==id
                o.flux = parseFloat(flux[o.text]);
            }
            return o;
        });
        return data;
    }

    function parse_flux_2(data, flux2) {
        data.reaction_paths = data.reaction_paths.map( function(o) {
            if (o.id in flux2 && o.flux) {
                o.flux = (parseFloat(flux2[o.id]) - o.flux);
            }
            return o;
        });
        data.reaction_labels = data.reaction_labels.map( function(o) {
            if (o.flux) o.flux1 = o.flux;
            else o.flux1 = 0;
            if (o.text in flux2) o.flux2 = parseFloat(flux2[o.text]);
            else o.flux2 = 0;
            o.flux = (o.flux2 - o.flux1);
            return o;
        });
        return data;
    }

    function parse_metabolites_1(data, metabolites) {
        skip_these_metabolites = []; //
        do_not_size_these_metabolites = ['nad','nadp','nadh','nadph','atp','adp','coa','accoa'];
        data.metabolite_circles = data.metabolite_circles.map( function(o) {
            if (o.id in metabolites && skip_these_metabolites.indexOf(o.id)==-1) {
                o.metabolite_concentration = parseFloat(metabolites[o.id])
                if (do_not_size_these_metabolites.indexOf(o.id)>=0) {
                    o.should_size = false;
                    o.should_color = true;
                } else {
                    o.should_size = true;
                    o.should_color = false;
                }
            }
            return o;
        });
        data.metabolite_labels = data.metabolite_labels.map( function(o) {
            if (o.text in metabolites) {
                o.metabolite_concentration = parseFloat(metabolites[o.text])
            }
            return o;
        });
        return data;
    }

    function parse_metabolites_2(data, metabolites) {
        data.metabolite_circles = data.metabolite_circles.map( function(o) {
            if (o.id in metabolites) {
                o.metabolite_deviation = parseFloat(metabolites[o.id])
            }
            return o;
        });
        data.metabolite_labels = data.metabolite_labels.map( function(o) {
            if (o.text in metabolites) {
                o.metabolite_deviation = parseFloat(metabolites[o.text])
            }
            return o;
        });
        return data;
    }

    function append_deviation_arcs(metabolite_circle_data) {
        arc_data = metabolite_circle_data.filter( function(o) {
            return (o.hasOwnProperty('metabolite_deviation') &&
                    o.hasOwnProperty('metabolite_concentration'));
        });
        var arc = d3.svg.arc()
            .startAngle(function(d) { return -d.metabolite_deviation/100/2*2*Math.PI; })
            .endAngle(function(d) { return d.metabolite_deviation/100/2*2*Math.PI; })
            .innerRadius(function(d) { return 0; })
            .outerRadius(function(d) {
                var s;
                if (d.should_size) s = maps.scale.size(maps.scale.metabolite_concentration(d.metabolite_concentration));
                else s = maps.scale.size(0);
                return s;
            });
        maps.svg.append("g")
            .attr("id", "metabolite-deviation-arcs")
            .selectAll("path")
            .data(arc_data)
            .enter().append("path")
            .attr('d', arc)
            .attr('style', "fill:black;stroke:none;opacity:0.4;")
            .attr("transform", function(d) {
                return "translate("+maps.scale.x(d.cx)+","+maps.scale.y(d.cy)+")";
            });
    }

    function append_reaction_paths(reaction_path_data) {
        maps.svg.append("g")
            .attr("id", "reaction-paths")
            .selectAll("path")
            .data(reaction_path_data)
            .enter().append("path")
            .attr("d", function(d) { return scale_path(d.d); })
            .attr("class", function(d) { return d.class })
            .attr("style", function(d) {
                var s = "", sc = maps.scale.flux;
                // .fill-arrow is for simpheny maps where the path surrounds line and
                // arrowhead
                // .line-arrow is for bigg maps were the line is a path and the
                // arrowhead is a marker
                if (d.class=="fill-arrow") sc = maps.scale.flux_fill;
                if (d.flux) {
                    s += "stroke-width:"+String(maps.scale.size(sc(Math.abs(d.flux))))+";";
                    s += "stroke:"+maps.scale.flux_color(Math.abs(d.flux))+";";
                    if (d.class=="fill-arrow") { s += "fill:"+maps.scale.flux_color(Math.abs(d.flux))+";"; }
                    else if (d.class=="line-arrow") { make_arrowhead_for_fill(); }
                    else s += "fill:none";
                }
                else if (maps.has_flux) {
                    s += "stroke-width:"+String(maps.scale.size(sc(0)))+";";
                    s += "stroke:"+maps.scale.flux_color(Math.abs(0))+";";
                    if (d.class=="fill-arrow") s += "fill:"+maps.scale.flux_color(0)+";";
                    else s += "fill:none";
                }
                else {
                    s += "stroke-width:"+String(maps.scale.size(1))+";";
                    s += "stroke:"+maps.defaults.path_color+";";
                    if (d.class=="fill-arrow") s += "fill:"+maps.defaults.path_color+";";
                    else s += "fill:none";
                }
                return s;
            })
            .style('marker-end', function (d) {
                if (!/end/.test(d.class)) return '';
                if (d.flux) return make_arrowhead_for_fill(maps.scale.flux_color(d.flux));
                else if (maps.has_flux) return make_arrowhead_for_fill(maps.scale.flux_color(0));
                else return "url(#end-triangle-path-color)";
            })
            .style('marker-start', function (d) {
                if (!/start/.test(d.class)) return '';
                if (d.flux) return make_arrowhead_for_fill(maps.scale.flux_color(d.flux));
                else if (maps.has_flux) return make_arrowhead_for_fill(maps.scale.flux_color(0));
                else return "url(#start-triangle-path-color)";
            });
    }

    function make_arrowhead_for_fill(fill) {
        d3.select('#markers').selectAll("marker"); //
        return ""
    }

    function scale_decimals(path, scale_fn, precision) {
        var str = d3.format("."+String(precision)+"f")
        path = path.replace(/([0-9.]+)/g, function (match, p1) {
            return str(scale_fn(parseFloat(p1)));
        });
        return path
    }

    function scale_path(path) {
        var x_fn = maps.scale.x; var y_fn = maps.scale.y;
        // TODO: scale arrow width
        var str = d3.format(".2f")
        path = path.replace(/(M|L)([0-9-.]+),?\s*([0-9-.]+)/g, function (match, p0, p1, p2) {
            return p0 + [str(x_fn(parseFloat(p1))), str(y_fn(parseFloat(p2)))].join(', ');
        });
        path = path.replace(/C([0-9-.]+),?\s*([0-9-.]+)\s*([0-9-.]+),?\s*([0-9-.]+)\s*([0-9-.]+),?\s*([0-9-.]+)/g, function (match, p1, p2, p3, p4, p5, p6) {
            return 'C'+str(x_fn(parseFloat(p1)))+','+str(y_fn(parseFloat(p2)))+' '+str(x_fn(parseFloat(p3)))+','+str(y_fn(parseFloat(p4)))+' '+ [str(x_fn(parseFloat(p5)))+','+str(y_fn(parseFloat(p6)))];
        });
        return path
    }

    return maps;
}(d3));

function download_map(div_name) {
    svg = d3.select("#" + div_name)[0][0].getElementsByTagName("svg")[0]
    var a = document.createElement('a'), xml, ev;
    a.download = 'map.svg'; // file name
    xml = (new XMLSerializer()).serializeToString(svg); // convert node to xml string
    a.setAttribute("href-lang", "image/svg+xml")
    a.href = 'data:image/svg+xml;base64,' + btoa(xml); // create data uri
    // <a> constructed, simulate mouse click on it
    ev = document.createEvent("MouseEvents");
    ev.initMouseEvent("click", true, false, self, 0, 0, 0, 0, 0, false, false, false, false, 0, null);
    a.dispatchEvent(ev);
}
