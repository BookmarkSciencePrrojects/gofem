#!/usr/bin/python

# Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

import subprocess

def Cmd(command, verbose=False, debug=False):
    if debug:
        print '=================================================='
        print cmd
        print '=================================================='
    spr = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = spr.stdout.read()
    err = spr.stderr.read().strip()
    if verbose:
        print out
        print err
    return out, err

pkgs = [
    ("ana",     "analytical solutions for comparisons"),
    ("shp",     "shape structures and quadrature points"),
    ("mdl/sld", "models for solids"),
    ("mdl/cnd", "models for liquid/gas conductivity in porous media"),
    ("mdl/lrm", "models for liquid retention in porous media"),
    ("mdl/fld", "models for fluids"),
    ("mdl/por", "models for porous media"),
    ("inp",     "input data structures. simulation, materials, meshes"),
    ("fem",     "finite element method"),
    ("out",     "results analyses and plotting"),
]

odir  = 'doc/'
idxfn = odir+'index.html'
licen = open('LICENSE', 'r').read()

def header(title):
    return """<html>
<head>
<meta http-equiv=\\"Content-Type\\" content=\\"text/html; charset=utf-8\\">
<title>%s</title>
<link type=\\"text/css\\" rel=\\"stylesheet\\" href=\\"static/style.css\\">
<script type=\\"text/javascript\\" src=\\"static/godocs.js\\"></script>
<style type=\\"text/css\\"></style>
</head>
<body>
<div id=\\"page\\">""" % title

def footer():
    return """</div><!-- page -->
<div id=\\"footer\\">
<br /><br />
<hr>
<pre class=\\"copyright\\">
%s</pre><!-- copyright -->
</div><!-- footer -->
</body>
</html>""" % licen

def pkgheader(pkg):
    return header('Gofem &ndash; package '+pkg[0]) + '<h1>Gofem &ndash; <b>%s</b> &ndash; %s</h1>' % (pkg[0], pkg[1])

def pkgitem(pkg):
    fnk = pkg[0].replace("/","-")
    return '<dd><a href=\\"xx%s.html\\"><b>%s</b>: %s</a></dd>' % (fnk, pkg[0], pkg[1])

Cmd('echo "'+header('Gofem &ndash; Documentation')+'" > '+idxfn)
Cmd('echo "<h1>Gofem &ndash; Documentation</h1>" >> '+idxfn)
Cmd('echo "<h2 id=\\"pkg-index\\">Index</h2>\n<div id=\\"manual-nav\\">\n<dl>" >> '+idxfn)

for pkg in pkgs:
    fnk = pkg[0].replace("/","-")
    fn = odir+'xx'+fnk+'.html'
    print fn
    Cmd('echo "'+pkgheader(pkg)+'" > '+fn)
    Cmd('godoc -html github.com/cpmech/gofem/'+pkg[0]+' >> '+fn)
    Cmd('echo "'+footer()+'" >> '+fn)
    Cmd('echo "'+pkgitem(pkg)+'" >> '+idxfn)

    # fix links
    Cmd("sed -i -e 's@/src/target@https://github.com/cpmech/gofem/blob/master/"+pkg[0]+"@g' "+fn+"")

Cmd('echo "</dl>\n</div><!-- manual-nav -->" >> '+idxfn)
Cmd('echo "'+footer()+'" >> '+idxfn)
