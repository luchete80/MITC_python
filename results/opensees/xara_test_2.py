import xara

model = xara.Model(ndm=2, ndf=2)

model.node(1, (0.0, 0.0))
model.node(2, (1.0, 0.0))
model.node(3, 0.0, 1.0,0)
model.material("Elastic", 1, 29e3, 0.3)
#model.element('ASDShellT3',1,  1,2,3,  1, '-corotational')
model.element('ASDShellT3',1,  1,2,3,  1, '--linear')
model.analysis("Static")
model.analyze(1)
