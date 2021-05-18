%% Create the HTML documentation files for the various classes & functions
addpath(genpath("../"))

%%  THzWaveformParameters
html = help2html("THzWaveformParameters",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen('THzWaveformParameters.html','w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% THzAntennaArray
html = help2html("THzAntennaArray",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen('THzAntennaArray.html','w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% THzScanner
html = help2html("THzScanner",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen('THzScanner.html','w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% THzTarget
html = help2html("THzTarget",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen('THzTarget.html','w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% THzImageReconstruction
html = help2html("THzImageReconstruction",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("THzImageReconstruction" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% nonuniform_theta_CSAR_XZ_BPA
html = help2html("nonuniform_theta_CSAR_XZ_BPA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("nonuniform_theta_CSAR_XZ_BPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% nonuniform_thetaY_CSAR_XYZ_BPA
html = help2html("nonuniform_thetaY_CSAR_XYZ_BPA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("nonuniform_thetaY_CSAR_XYZ_BPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% nonuniform_XY_SAR_XY_BPA
html = help2html("nonuniform_XY_SAR_XY_BPA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("nonuniform_XY_SAR_XY_BPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% nonuniform_XY_SAR_XYZ_BPA
html = help2html("nonuniform_XY_SAR_XYZ_BPA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("nonuniform_XY_SAR_XYZ_BPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% nonuniform_Y_SAR_YZ_BPA
html = help2html("nonuniform_Y_SAR_YZ_BPA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("nonuniform_Y_SAR_YZ_BPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% reconstructionAlgorithmTemplate
html = help2html("reconstructionAlgorithmTemplate",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("reconstructionAlgorithmTemplate" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_theta_CSAR_XZ_PFA
html = help2html("uniform_theta_CSAR_XZ_PFA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("uniform_theta_CSAR_XZ_PFA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_thetaY_CSAR_XYZ_PFA
html = help2html("uniform_thetaY_CSAR_XYZ_PFA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("uniform_thetaY_CSAR_XYZ_PFA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_XY_SAR_XY_FFT
html = help2html("uniform_XY_SAR_XY_FFT",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("uniform_XY_SAR_XY_FFT" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_XY_SAR_XYZ_RMA
html = help2html("uniform_XY_SAR_XYZ_RMA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("uniform_XY_SAR_XYZ_RMA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_Y_SAR_Y_FFT
html = help2html("uniform_Y_SAR_Y_FFT",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("uniform_Y_SAR_Y_FFT" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% uniform_Y_SAR_YZ_RMA
html = help2html("uniform_Y_SAR_YZ_RMA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("uniform_Y_SAR_YZ_RMA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% computeIdeal2D
html = help2html("computeIdeal2D",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("computeIdeal2D" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% computeIdeal3D
html = help2html("computeIdeal3D",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("computeIdeal3D" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% displayImage1D
html = help2html("displayImage1D",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("displayImage1D" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% displayImage2D
html = help2html("displayImage2D",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("displayImage2D" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% displayImage3D
html = help2html("displayImage3D",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("displayImage3D" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% fastBPA
html = help2html("fastBPA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("fastBPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% getFields
html = help2html("getFields",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("getFields" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% isApp
html = help2html("isApp",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("isApp" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% isAppEPC
html = help2html("isAppEPC",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("isAppEPC" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% mediumBPA
html = help2html("mediumBPA",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("mediumBPA" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% mult2mono
html = help2html("mult2mono",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("mult2mono" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% plotXYdB
html = help2html("plotXYdB",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("plotXYdB" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% plotXYZdB
html = help2html("plotXYZdB",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("plotXYZdB" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);

%% showErrorMessage
html = help2html("showErrorMessage",[],"-doc");
html = strrep(html,'https://www.mathworks.com/help/releases/R2021a/includes/product/css/helpwin.css','helpwin2.css');
fid = fopen("showErrorMessage" + ".html",'w');
fprintf(fid,'%s',string(html));
fclose(fid);
