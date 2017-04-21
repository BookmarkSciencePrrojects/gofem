#!/usr/bin/python

# Copyright 2015 The Gosl Authors. All rights reserved.
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
    ("ana"              , "analytical solutions"),
    ("shp"              , "shape (interpolation) structures and quadrature points"),
    ("mdl/generic"      , "generic models (placeholder for parameters set)"),
    ("mdl/solid"        , "models for solids"),
    ("mdl/fluid"        , "models for fluids (liquid / gas)"),
    ("mdl/conduct"      , "models for liquid conductivity within porous media"),
    ("mdl/retention"    , "models for liquid retention within porous media"),
    ("mdl/diffusion"    , "models for diffusion applications"),
    ("mdl/thermomech"   , "models for thermo-mechanical applications"),
    ("mdl/porous"       , "models for porous media (TPM-based)"),
    ("inp"              , "input data (.sim = simulation, .mat = materials, .msh = meshes)"),
    ("ele"              , "finite elements"),
    ("ele/solid"        , "elements for solid mechanics"),
    ("ele/seepage"      , "elements for seepage problems (with liquid and/or gases)"),
    ("ele/diffusion"    , "elements for diffusion(-like) problems"),
    ("ele/thermomech"   , "elements for thermo-mechanical applications"),
    ("ele/porous"       , "elements for porous media simulations (TPM)"),
    ("fem"              , "finite element method (main, domain, solver, ...)"),
    ("tests"            , "(unit) tests of complete simulations"),
    ("tests/solid"      , "tests of solid mechanics applications"),
    ("tests/seepage"    , "tests of seepage problems"),
    ("tests/diffusion"  , "tests of diffusion problems"),
    ("tests/thermomech" , "tests of thermo-mechanical applications"),
    ("tests/porous"     , "tests of porous media simulations"),
    ("out"              , "output routines (post-processing and plotting)"),
]

odir  = 'doc/'
idxfn = odir+'index.html'
licen = open('LICENSE', 'r').read()

def header(title):
    return """<html>
<head>
<meta http-equiv=\\"Content-Type\\" content=\\"text/html; charset=utf-8\\">
<meta name=\\"viewport\\" content=\\"width=device-width, initial-scale=1\\">
<meta name=\\"theme-color\\" content=\\"#375EAB\\">
<title>%s</title>
<link type=\\"text/css\\" rel=\\"stylesheet\\" href=\\"static/style.css\\">
<script type=\\"text/javascript\\" src=\\"static/godocs.js\\"></script>
<style type=\\"text/css\\"></style>
</head>
<body>
<div id=\\"page\\" class=\\wide\\">
<div class=\\"container\\">
""" % title

def footer():
    return """
<div id=\\"footer\\">
<br /><br />
<hr>
<pre class=\\"copyright\\">
%s</pre><!-- copyright -->
</div><!-- footer -->

</div><!-- container -->
</div><!-- page -->
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
    Cmd("sed -i -e 's@/src/github.com/cpmech/gofem/@https://github.com/cpmech/gofem/blob/master/@g' "+fn+"")

    # fix links to subdirectories (harder to automate)
    subdirs = []
    if pkg[0] == "mdl/solid":
        subdirs = ["data"]

    if pkg[0] == "mdl/porous":
        subdirs = ["data"]

    if pkg[0] == "inp":
        subdirs = ["data"]

    if pkg[0] == "ele":
        subdirs = ["diffusion", "porous", "seepage", "solid", "thermomech"]

    if pkg[0] == "ele/solid":
        subdirs = ["data"]

    if pkg[0] == "ele/seepage":
        subdirs = ["data"]

    if pkg[0] == "ele/diffusion":
        subdirs = ["data"]

    if pkg[0] == "ele/porous":
        subdirs = ["data"]

    if pkg[0] == "fem":
        subdirs = ["data"]

    if pkg[0] == "tests":
        subdirs = ["data", "diffusion", "porous", "seepage", "solid", "thermomech"]

    if pkg[0] == "tests/solid":
        subdirs = ["cmp", "data"]

    if pkg[0] == "tests/seepage":
        subdirs = ["data"]

    if pkg[0] == "tests/diffusion":
        subdirs = ["data"]

    if pkg[0] == "tests/thermomech":
        subdirs = ["data"]

    if pkg[0] == "tests/porous":
        subdirs = ["data"]

    if pkg[0] == "out":
        subdirs = ["cmp", "data"]

    for subdir in subdirs:
        Cmd("sed -i -e 's@<a href=\""+subdir+"/\">@<a href=\"https://github.com/cpmech/gofem/tree/master/"+pkg[0]+"/"+subdir+"\">@g' "+fn)


    # remove link to parent directory
    Cmd("sed -i -e 's@<td colspan=\"2\"><a href=\"..\">..</a></td>@<td></td>@g' "+fn)

Cmd('echo "</dl>\n</div><!-- manual-nav -->" >> '+idxfn)
Cmd('echo "'+footer()+'" >> '+idxfn)
