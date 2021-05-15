%% Create the HTML documentation files for the various classes & functions
addpath(genpath("../"))

%%  fmcwChirpParameters
html = help2html("fmcwChirpParameters",[],"-doc");
fid = fopen('fmcwChirpParameters.html','w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% sarAntennaArray
html = help2html("sarAntennaArray",[],"-doc");
fid = fopen('sarAntennaArray.html','w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% sarScenario
html = help2html("sarScenario",[],"-doc");
fid = fopen('sarScenario.html','w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% sarTarget
html = help2html("sarTarget",[],"-doc");
fid = fopen('sarTarget.html','w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% sarImage
html = help2html("sarImage",[],"-doc");
fid = fopen("sarImage" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% nonuniform_theta_CSAR_XZ_BPA
html = help2html("nonuniform_theta_CSAR_XZ_BPA",[],"-doc");
fid = fopen("nonuniform_theta_CSAR_XZ_BPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% nonuniform_thetaY_CSAR_XYZ_BPA
html = help2html("nonuniform_thetaY_CSAR_XYZ_BPA",[],"-doc");
fid = fopen("nonuniform_thetaY_CSAR_XYZ_BPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% nonuniform_XY_SAR_XY_BPA
html = help2html("nonuniform_XY_SAR_XY_BPA",[],"-doc");
fid = fopen("nonuniform_XY_SAR_XY_BPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% nonuniform_XY_SAR_XYZ_BPA
html = help2html("nonuniform_XY_SAR_XYZ_BPA",[],"-doc");
fid = fopen("nonuniform_XY_SAR_XYZ_BPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% reconstructionAlgorithmTemplate
html = help2html("reconstructionAlgorithmTemplate",[],"-doc");
fid = fopen("reconstructionAlgorithmTemplate" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_theta_CSAR_XZ_PFA
html = help2html("uniform_theta_CSAR_XZ_PFA",[],"-doc");
fid = fopen("uniform_theta_CSAR_XZ_PFA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_thetaY_CSAR_XYZ_PFA
html = help2html("uniform_thetaY_CSAR_XYZ_PFA",[],"-doc");
fid = fopen("uniform_thetaY_CSAR_XYZ_PFA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_XY_SAR_XY_FFT
html = help2html("uniform_XY_SAR_XY_FFT",[],"-doc");
fid = fopen("uniform_XY_SAR_XY_FFT" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_XY_SAR_XYZ_RMA
html = help2html("uniform_XY_SAR_XYZ_RMA",[],"-doc");
fid = fopen("uniform_XY_SAR_XYZ_RMA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_Y_SAR_Y_FFT
html = help2html("uniform_Y_SAR_Y_FFT",[],"-doc");
fid = fopen("uniform_Y_SAR_Y_FFT" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_Y_SAR_YZ_RMA
html = help2html("uniform_Y_SAR_YZ_RMA",[],"-doc");
fid = fopen("uniform_Y_SAR_YZ_RMA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% displayImage1D
html = help2html("displayImage1D",[],"-doc");
fid = fopen("displayImage1D" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% displayImage2D
html = help2html("displayImage2D",[],"-doc");
fid = fopen("displayImage2D" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% displayImage3D
html = help2html("displayImage3D",[],"-doc");
fid = fopen("displayImage3D" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% fastBPA
html = help2html("fastBPA",[],"-doc");
fid = fopen("fastBPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% mediumBPA
html = help2html("mediumBPA",[],"-doc");
fid = fopen("mediumBPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% mult2mono
html = help2html("mult2mono",[],"-doc");
fid = fopen("mult2mono" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% plotXYdB
html = help2html("plotXYdB",[],"-doc");
fid = fopen("plotXYdB" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% plotXYZdB
html = help2html("plotXYZdB",[],"-doc");
fid = fopen("plotXYZdB" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);
