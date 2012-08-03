vtkRenderWindow renWin
renWin SetSize 600 300
renWin SetPosition 600 400

vtkRenderer ren1
ren1 SetViewport 0.0 0.0 0.5 1.0
ren1 SetBackground 0.8 0.4 0.2
renWin AddRenderer ren1

vtkRenderer ren2
ren2 SetViewport 0.5 0.0 1.0 1.0
ren2 SetBackground 0.1 0.2 0.4
renWin AddRenderer ren2

renWin Render

vtkCubeSource cubeData
vtkPolyDataMapper cubeMapper
cubeMapper SetInput [cubeData GetOutput]

vtkActor cubeActor
cubeActor SetMapper cubeMapper

ren1 AddViewProp cubeActor
ren1 ResetCamera
renWin Render

cubeActor RotateX 30.0
cubeActor RotateY 20.0
[cubeActor GetProperty] SetColor 1.0 0.7 0.7

renWin Render

vtkSLCReader negReader
negReader SetFileName "neghip.slc"

vtkVolumeTextureMapper2D negMapper
negMapper SetInput [negReader GetOutput]

vtkPiecewiseFunction negOpacity
negOpacity AddPoint 0 0.0
negOpacity AddPoint 255 0.2

vtkColorTransferFunction negColor
negColor AddRGBPoint  64 1.0 0.0 0.0
negColor AddRGBPoint 128 0.0 0.0 1.0
negColor AddRGBPoint 196 0.0 1.0 0.0

vtkVolumeProperty negProperty
negProperty SetColor negColor
negProperty SetScalarOpacity negOpacity

vtkVolume negVolume
negVolume SetMapper negMapper
negVolume SetProperty negProperty

ren2 AddViewProp negVolume
ren2 ResetCamera
renWin Render

#vtkPolyDataReader posReader
#posReader SetFileName "poshipsurface.vtk"

#vtkPolyDataMapper posMapper
#posMapper SetInput [posReader GetOutput]

#vtkActor posActor
#posActor SetMapper posMapper

#ren2 AddViewProp posActor
#renWin Render

vtkTextMapper titleMapper
titleMapper SetInput "This is a pink cube"
#titleMapper SetJustificationToCenter

vtkActor2D titleActor
titleActor SetMapper titleMapper
[titleActor GetProperty] SetColor 1 1 0
set pc [titleActor GetPositionCoordinate]
$pc SetCoordinateSystemToNormalizedViewport
$pc SetValue 0.5 0.92

ren1 AddViewProp titleActor
renWin Render

proc changeTitle {} {
	titleMapper SetInput [.top.entry get]
	renWin Render
}	

toplevel .top
entry .top.entry
.top.entry insert 0 "Type Title for left Image here"
pack .top.entry
bind .top.entry <Return> changeTitle