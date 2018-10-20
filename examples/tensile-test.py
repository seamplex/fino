# state file generated using paraview version 5.1.2

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1215, 822]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 0
renderView1.CameraPosition = [119.85527680930403, 83.5687295595565, 215.494572652346]
renderView1.CameraFocalPoint = [25.0983059214126, -7.58028478306764, -7.56438553653127]
renderView1.CameraViewUp = [-0.09535621980896813, 0.9350079031724332, -0.34156611709716406]
renderView1.CameraParallelScale = 81.0879152525208
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.Background2 = [0.0, 0.0, 0.164705882352941]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
tensiletestvtk = LegacyVTKReader(FileNames=['/home/gtheler/codigos/wasora-suite-0.5/fino/examples/tensile-test.vtk'])

# create a new 'Warp By Vector'
warpByVector1 = WarpByVector(Input=tensiletestvtk)
warpByVector1.Vectors = ['POINTS', 'u-v-w']
warpByVector1.ScaleFactor = 50.0

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get color transfer function/color map for 'sigma'
sigmaLUT = GetColorTransferFunction('sigma')
sigmaLUT.LockDataRange = 1
sigmaLUT.RGBPoints = [0.0, 0.0, 0.0, 0.5625, 66.6666, 0.0, 0.0, 1.0, 219.0477, 0.0, 1.0, 1.0, 295.2381, 0.5, 1.0, 0.5, 371.4285, 1.0, 1.0, 0.0, 523.8095999999999, 1.0, 0.0, 0.0, 600.0, 0.5, 0.0, 0.0]
sigmaLUT.ColorSpace = 'RGB'
sigmaLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
sigmaLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'sigma'
sigmaPWF = GetOpacityTransferFunction('sigma')
sigmaPWF.Points = [0.0, 0.0, 0.5, 0.0, 600.0, 1.0, 0.5, 0.0]
sigmaPWF.ScalarRangeInitialized = 1

# get color transfer function/color map for 'u-v-w'
uvwLUT = GetColorTransferFunction('u-v-w')
uvwLUT.RGBPoints = [0.0, 0.23, 0.299, 0.754, 0.150814646577962, 0.865, 0.865, 0.865, 0.301629293155924, 0.706, 0.016, 0.15]
uvwLUT.NanColor = [0.25, 0.0, 0.0]
uvwLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'u-v-w'
uvwPWF = GetOpacityTransferFunction('u-v-w')
uvwPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.301629293155924, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from tensiletestvtk
tensiletestvtkDisplay = Show(tensiletestvtk, renderView1)
# trace defaults for the display properties.
tensiletestvtkDisplay.ColorArrayName = ['POINTS', 'sigma']
tensiletestvtkDisplay.LookupTable = sigmaLUT
tensiletestvtkDisplay.Specular = 0.1
tensiletestvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
tensiletestvtkDisplay.GlyphType = 'Arrow'
tensiletestvtkDisplay.SelectionPointLabelColor = [0.5, 0.5, 0.5]
tensiletestvtkDisplay.ScalarOpacityUnitDistance = 10.910678224057

# show color legend
tensiletestvtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data from warpByVector1
warpByVector1Display = Show(warpByVector1, renderView1)
# trace defaults for the display properties.
warpByVector1Display.Representation = 'Wireframe'
warpByVector1Display.AmbientColor = [0.0, 0.0, 0.0]
warpByVector1Display.ColorArrayName = ['POINTS', '']
warpByVector1Display.LookupTable = uvwLUT
warpByVector1Display.Specular = 0.1
warpByVector1Display.OSPRayScaleFunction = 'PiecewiseFunction'
warpByVector1Display.GlyphType = 'Arrow'
warpByVector1Display.SelectionPointLabelColor = [0.5, 0.5, 0.5]
warpByVector1Display.ScalarOpacityUnitDistance = 10.930702054107

# setup the color legend parameters for each legend in this view

# get color legend/bar for sigmaLUT in view renderView1
sigmaLUTColorBar = GetScalarBar(sigmaLUT, renderView1)
sigmaLUTColorBar.Position = [0.87, 0.25]
sigmaLUTColorBar.Position2 = [0.13, 0.5]
sigmaLUTColorBar.Title = 'Von Mises'
sigmaLUTColorBar.ComponentTitle = '[MPa]'
sigmaLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
sigmaLUTColorBar.TitleFontFamily = 'Times'
sigmaLUTColorBar.TitleFontSize = 8
sigmaLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
sigmaLUTColorBar.LabelFontSize = 6
sigmaLUTColorBar.AutomaticLabelFormat = 0
sigmaLUTColorBar.LabelFormat = '%g'
sigmaLUTColorBar.RangeLabelFormat = '%g'

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(tensiletestvtk)
# ----------------------------------------------------------------