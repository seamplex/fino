try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

a1_T_PiecewiseFunction = CreatePiecewiseFunction( Points=[397.319, 0.0, 0.5, 0.0, 593.392, 1.0, 0.5, 0.0] )

a1_T_PVLookupTable = GetLookupTableForArray( "T", 1, RGBPoints=[397.319, 0.0, 0.0, 1.0, 593.392, 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.498039, 0.498039, 0.498039], ScalarOpacityFunction=a1_T_PiecewiseFunction, ColorSpace='HSV', ScalarRangeInitialized=1.0 )

RenderView1 = CreateRenderView()
RenderView1.CameraViewUp = [0.35476473200364605, 0.8890321694918081, -0.28941974109426305]
RenderView1.CacheKey = 0.0
RenderView1.UseGradientBackground = 1
RenderView1.StereoType = 0
RenderView1.UseLight = 1
RenderView1.StereoRender = 0
RenderView1.CameraPosition = [-113.971945306639, 54.5740603002526, 109.656956416655]
RenderView1.LightSwitch = 0
RenderView1.OrientationAxesVisibility = 0
RenderView1.Background2 = [0.9333333333333333, 0.9333333333333333, 0.9333333333333333]
RenderView1.CameraClippingRange = [99.86536979083834, 307.35797957811883]
RenderView1.StereoCapableWindow = 0
RenderView1.Background = [0.8, 0.8, 0.8]
RenderView1.CameraFocalPoint = [29.6766491130484, -45.1826789668781, -20.6917719992069]
RenderView1.CenterAxesVisibility = 0
RenderView1.CameraParallelScale = 56.4540742550969
RenderView1.CenterOfRotation = [0.0, -33.75, 0.0]

ScalarBarWidgetRepresentation1 = CreateScalarBar( LabelFormat='%-#6.2f', Title='T', Enabled=1, LabelFontSize=10, LabelColor=[0.0, 0.0, 0.5019607843137255], LabelFontFamily='Times', LookupTable=a1_T_PVLookupTable, TitleFontSize=10, TitleColor=[0.0, 0.0, 0.5019607843137255], TitleFontFamily='Times' )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)

enginepiston_vtk = LegacyVTKReader( guiName="engine-piston.vtk", FileNames=['/home/gtheler/codigos/wasora-suite/fino/examples/engine-piston.vtk'] )

Clip1 = Clip( guiName="Clip1", Scalars=['POINTS', 'T'], ClipType="Plane", Value=495.3515 )
Clip1.ClipType.Origin = [0.0, -33.75, 0.0]

SetActiveSource(enginepiston_vtk)
DataRepresentation1 = Show()
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation1.SelectionPointFieldDataArrayName = 'T'
DataRepresentation1.ScalarOpacityFunction = a1_T_PiecewiseFunction
DataRepresentation1.ColorArrayName = ('POINT_DATA', 'T')
DataRepresentation1.ScalarOpacityUnitDistance = 3.24187864957382
DataRepresentation1.Visibility = 0
DataRepresentation1.LookupTable = a1_T_PVLookupTable
DataRepresentation1.ScaleFactor = 6.75

SetActiveSource(Clip1)
DataRepresentation2 = Show()
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.SelectionPointFieldDataArrayName = 'T'
DataRepresentation2.ScalarOpacityFunction = a1_T_PiecewiseFunction
DataRepresentation2.ColorArrayName = ('POINT_DATA', 'T')
DataRepresentation2.ScalarOpacityUnitDistance = 3.52468522839747
DataRepresentation2.LookupTable = a1_T_PVLookupTable
DataRepresentation2.ScaleFactor = 6.75


GetRenderView().ViewSize = [1200,800]
Render()

WriteImage("engine-piston.png")
