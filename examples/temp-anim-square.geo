Merge "temp-square.msh";

end_time = 2000;
dt = 10;
For time In {0:end_time:dt}
  View[0].Time = time;
  Print Sprintf("temp-square-rough-%03g.png", time/dt);
//  Sleep 0.01;
  Draw;
EndFor

General.Terminal = 1;
Printf("all frames dumped, now run");
Printf("ffmpeg -y -framerate 20 -f image2 -i temp-square-rough-%%03d.png temp-square-rough.mp4");
Printf("to get a video");
Exit;
