
SSPipeline
==========

SSPipeline stands for [Slicer](http://www.slicer.org) Slice Pipeline.

It is a small experiment that should help understanding how images are rendered in the Slice slice views.

It implements the pipeline components required to render the scalar volume MRHead as a Axial foreground layer.

Prerequisites
-------------

* Slicer build tree. See [here](http://wiki.slicer.org/slicerWiki/index.php/Documentation/Nightly/Developers/Build_Instructions) for more details.


Building
--------

```
git clone git://github.com/jcfr/SSPipeline

mkdir SSPipeline-build && cd SSPipeline-build

cmake -DCMAKE_BUILD_TYPE:STRING=Debug \
  -DSlicer_DIR:PATH=/home/jchris/Projects/Slicer-SuperBuild-Debug-Qt485/Slicer-build/ \
  ../SSPipeline
```

References
----------

* http://wiki.slicer.org/slicerWiki/index.php/Documentation/Nightly/Developers/MRML#Slice_view_pipeline
