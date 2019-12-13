# state file generated using paraview version 5.4.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1205, 807]
renderView1.AnnotationColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.StereoType = 0
renderView1.CameraPosition = [142.80065864526452, 114.67039791632159, 185.69279953353038]
renderView1.CameraFocalPoint = [18.52862152159883, 16.47005370719967, -19.137565310482255]
renderView1.CameraViewUp = [-0.25016855087924833, 0.9234534904032745, -0.2909456086848301]
renderView1.CameraParallelScale = 81.08791525252083
renderView1.Background = [1.0, 1.0, 1.0]

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
tensiletestvtk = LegacyVTKReader(FileNames=['/home/gtheler/codigos/wasora-suite/fino/examples/tensile-test.vtk'])

# create a new 'Warp By Vector'
warpByVector1 = WarpByVector(Input=tensiletestvtk)
warpByVector1.Vectors = ['POINTS', 'u-v-w']
warpByVector1.ScaleFactor = 200.0

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get color transfer function/color map for 'sigma'
sigmaLUT = GetColorTransferFunction('sigma')
sigmaLUT.RGBPoints = [3.68617, 0.32549, 0.14902, 0.960784, 14.125201556785, 0.297047, 0.375586, 0.963836, 24.564311145485, 0.180302, 0.536818, 0.964627, 35.00334270227, 0.1302, 0.649207, 0.929647, 45.442374259055, 0.0445143, 0.749654, 0.855998, 55.88140581583999, 0.0271325, 0.830713, 0.721527, 66.32051540453999, 0.259504, 0.866145, 0.543555, 76.759531354942, 0.428364, 0.890725, 0.329819, 87.1985863213015, 0.568503, 0.898508, 0.187623, 97.63761007489501, 0.738259, 0.890317, 0.0825461, 108.076719663595, 0.84546, 0.86136, 0.0147555, 118.51575122038001, 0.912191, 0.808018, 0.0, 128.95478277716498, 0.962848, 0.710445, 0.0, 139.39381433395, 0.999469, 0.600258, 0.0176284, 149.83292392265, 0.994156, 0.445975, 0.193912, 159.75, 0.980407, 0.247105, 0.262699]
sigmaLUT.ColorSpace = 'Lab'
sigmaLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'sigma'
sigmaPWF = GetOpacityTransferFunction('sigma')
sigmaPWF.Points = [3.68617, 0.0, 0.5, 0.0, 159.75, 1.0, 0.5, 0.0]
sigmaPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from tensiletestvtk
tensiletestvtkDisplay = Show(tensiletestvtk, renderView1)
# trace defaults for the display properties.
tensiletestvtkDisplay.Representation = 'Wireframe'
tensiletestvtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
tensiletestvtkDisplay.ColorArrayName = ['POINTS', '']
tensiletestvtkDisplay.DiffuseColor = [1.0, 0.6666666666666666, 1.0]
tensiletestvtkDisplay.Opacity = 0.5
tensiletestvtkDisplay.BackfaceDiffuseColor = [1.0, 0.6666666666666666, 1.0]
tensiletestvtkDisplay.OSPRayScaleArray = 'sigma'
tensiletestvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
tensiletestvtkDisplay.SelectOrientationVectors = 'u-v-w'
tensiletestvtkDisplay.ScaleFactor = 16.0
tensiletestvtkDisplay.SelectScaleArray = 'sigma'
tensiletestvtkDisplay.GlyphType = 'Arrow'
tensiletestvtkDisplay.GlyphTableIndexArray = 'sigma'
tensiletestvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
tensiletestvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
tensiletestvtkDisplay.ScalarOpacityUnitDistance = 5.362711587252902
tensiletestvtkDisplay.GaussianRadius = 8.0
tensiletestvtkDisplay.SetScaleArray = ['POINTS', 'sigma']
tensiletestvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
tensiletestvtkDisplay.OpacityArray = ['POINTS', 'sigma']
tensiletestvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
tensiletestvtkDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
tensiletestvtkDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
tensiletestvtkDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
tensiletestvtkDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
tensiletestvtkDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
tensiletestvtkDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
tensiletestvtkDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
tensiletestvtkDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
tensiletestvtkDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
tensiletestvtkDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
tensiletestvtkDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show data from warpByVector1
warpByVector1Display = Show(warpByVector1, renderView1)
# trace defaults for the display properties.
warpByVector1Display.Representation = 'Surface'
warpByVector1Display.AmbientColor = [0.0, 0.0, 0.0]
warpByVector1Display.ColorArrayName = ['POINTS', 'sigma']
warpByVector1Display.DiffuseColor = [1.0, 0.6666666666666666, 1.0]
warpByVector1Display.LookupTable = sigmaLUT
warpByVector1Display.BackfaceDiffuseColor = [1.0, 0.6666666666666666, 1.0]
warpByVector1Display.OSPRayScaleArray = 'sigma'
warpByVector1Display.OSPRayScaleFunction = 'PiecewiseFunction'
warpByVector1Display.SelectOrientationVectors = 'u-v-w'
warpByVector1Display.ScaleFactor = 17.510758
warpByVector1Display.SelectScaleArray = 'sigma'
warpByVector1Display.GlyphType = 'Arrow'
warpByVector1Display.GlyphTableIndexArray = 'sigma'
warpByVector1Display.DataAxesGrid = 'GridAxesRepresentation'
warpByVector1Display.PolarAxes = 'PolarAxesRepresentation'
warpByVector1Display.ScalarOpacityFunction = sigmaPWF
warpByVector1Display.ScalarOpacityUnitDistance = 5.856143502147454
warpByVector1Display.GaussianRadius = 8.755379
warpByVector1Display.SetScaleArray = ['POINTS', 'sigma']
warpByVector1Display.ScaleTransferFunction = 'PiecewiseFunction'
warpByVector1Display.OpacityArray = ['POINTS', 'sigma']
warpByVector1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
warpByVector1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
warpByVector1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
warpByVector1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
warpByVector1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
warpByVector1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
warpByVector1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
warpByVector1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
warpByVector1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
warpByVector1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
warpByVector1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
warpByVector1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show color legend
warpByVector1Display.SetScalarBarVisibility(renderView1, True)

# setup the color legend parameters for each legend in this view

# get color legend/bar for sigmaLUT in view renderView1
sigmaLUTColorBar = GetScalarBar(sigmaLUT, renderView1)
sigmaLUTColorBar.Title = 'sigma'
sigmaLUTColorBar.ComponentTitle = ''
sigmaLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
sigmaLUTColorBar.LabelColor = [0.0, 0.0, 0.0]

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(tensiletestvtk)
# ----------------------------------------------------------------