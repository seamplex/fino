try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

a1_sigma_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 480.0, 1.0, 0.5, 0.0] )

a3_uvw_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 0.3016292931559241, 1.0, 0.5, 0.0] )

a1_sigma_PVLookupTable = GetLookupTableForArray( "sigma", 1, RGBPoints=[0.0, 0.0, 0.0, 1.0, 480.0, 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.498039, 0.498039, 0.498039], ScalarOpacityFunction=a1_sigma_PiecewiseFunction, ColorSpace='HSV', ScalarRangeInitialized=1.0, LockScalarRange=1 )

a3_uvw_PVLookupTable = GetLookupTableForArray( "u-v-w", 3, RGBPoints=[0.0, 0.23, 0.299, 0.754, 0.15081464657796206, 0.865, 0.865, 0.865, 0.3016292931559241, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ScalarOpacityFunction=a3_uvw_PiecewiseFunction, ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

RenderView1 = CreateRenderView()
RenderView1.CameraViewUp = [-0.09535621980896805, 0.9350079031724333, -0.34156611709716383]
RenderView1.CacheKey = 0.0
RenderView1.StereoType = 0
RenderView1.UseLight = 1
RenderView1.StereoRender = 0
RenderView1.CameraPosition = [119.85527680930373, 83.56872955955654, 215.49457265234577]
RenderView1.LightSwitch = 0
RenderView1.OrientationAxesVisibility = 0
RenderView1.Background2 = [0.0, 0.0, 0.16470588235294117]
RenderView1.CameraClippingRange = [175.80556336151966, 357.445137619535]
RenderView1.StereoCapableWindow = 0
RenderView1.Background = [1.0, 1.0, 1.0]
RenderView1.CameraFocalPoint = [25.098305921412578, -7.580284783067641, -7.564385536531268]
RenderView1.CenterAxesVisibility = 0
RenderView1.CameraParallelScale = 81.08791525252083

ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='sigma', Enabled=1, LabelFontSize=6, LabelColor=[0.0, 0.0, 0.0], LookupTable=a1_sigma_PVLookupTable, TitleFontSize=8, TitleColor=[0.0, 0.0, 0.0], TitleFontFamily='Times' )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)

tensiletest_vtk = LegacyVTKReader( guiName="tensile-test.vtk", FileNames=['tensile-test.vtk'] )

WarpByVector1 = WarpByVector( guiName="WarpByVector1", Vectors=['POINTS', 'u-v-w'], ScaleFactor=50.0 )

SetActiveSource(tensiletest_vtk)
DataRepresentation1 = Show()
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation1.SelectionPointFieldDataArrayName = 'sigma'
DataRepresentation1.ScalarOpacityFunction = a1_sigma_PiecewiseFunction
DataRepresentation1.ColorArrayName = ('POINT_DATA', 'sigma')
DataRepresentation1.ScalarOpacityUnitDistance = 10.91067822405703
DataRepresentation1.LookupTable = a1_sigma_PVLookupTable
DataRepresentation1.ScaleFactor = 16.0

SetActiveSource(WarpByVector1)
DataRepresentation2 = Show()
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.SelectionPointFieldDataArrayName = 'sigma'
DataRepresentation2.ScalarOpacityFunction = a3_uvw_PiecewiseFunction
DataRepresentation2.ScalarOpacityUnitDistance = 10.930702054106991
DataRepresentation2.AmbientColor = [0.0, 0.0, 0.0]
DataRepresentation2.LookupTable = a3_uvw_PVLookupTable
DataRepresentation2.Representation = 'Wireframe'
DataRepresentation2.ScaleFactor = 16.0301352

Render()
