//#define SHOW_DEBUG_IMAGE
//# define MY_DEBUG
#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
//#include <stdio.h>
//#include <tchar.h>

#include <cstdio>	// Used for "printf"
#include <string>	// Used for C++ strings
#include <iostream>	// Used for C++ cout print statements
#include <fstream>
//#include <cmath>	// Used to calculate square-root for statistics

// Include OpenCV libraries
#include <cv.h>
#include <cvaux.h>
#include <cxcore.h>
#include <highgui.h>

using namespace std;

// Various color types for detected shirt colors.
enum                             {cBLACK=0,cWHITE, cGREY, cRED, cORANGE, cYELLOW, cGREEN, cAQUA, cBLUE, cPURPLE, cPINK,  NUM_COLOR_TYPES};
char* sCTypes[NUM_COLOR_TYPES] =    {"Black", "White","Grey","Red","Orange","Yellow","Green","Aqua","Blue","Purple","Pink"};
uchar cCTHue[NUM_COLOR_TYPES] =    {0,       0,        0,      0,      20,       30,       55,     85,    115,    138,  161};
uchar cCTSat[NUM_COLOR_TYPES] =     {0,       0,        0,      255,   255,     255,     255,   255,   255,    255,  255};
uchar cCTVal[NUM_COLOR_TYPES] =     {0,       255,    120,   255,   255,     255,     255,   255,   255,    255,  255};

int getPixelColorBlock(int H,int S,int V,int num_color_types){
	int color=-1;
    int color_block = num_color_types  - 3; // 0 black 1 white 2 gray
    int block_width = 180 / color_block;
	if (V < 75)
		color = 0;
	else if (V > 180 && S < 47)
		color = 1;
	else if (V < 175 && S < 63)
		color = 2;
	else {	// Is a color
        while(H>0){
            color++;
            H -= block_width;
        }
        color +=3;
    }
    return color;
}

int getPixelColorType(int H, int S, int V)
{
	int color;
	if (V < 75)
		color = cBLACK;
	else if (V > 180 && S < 47)
		color = cWHITE;
	else if (V < 175 && S < 63)
		color = cGREY;
	else {	// Is a color
		if (H < 14)
			color = cRED;
		else if (H < 25)
			color = cORANGE;
		else if (H < 34)
			color = cYELLOW;
		else if (H < 73)
			color = cGREEN;
		else if (H < 102)
			color = cAQUA;
		else if (H < 127)
			color = cBLUE;
		else if (H < 149)
			color = cPURPLE;
		else if (H < 165)
			color = cPINK;
		else	// full circle 
			color = cRED;	// back to Red
	}
	return color;
}

void getColorHistogram(IplImage* img, const CvRect& croprect, double colorHistogram[], int num_color_types, int type=0);

// hist1 -> center hist2 -> background
int combineHistogram(int hist_size, double hist1[], double hist2[]){
    int hist1_max=-1;
    int hist2_max=-1;
    int hist1_max_new = -1;
    double max_color_v  = -1;
   //
    max_color_v  = -1;
    for(int i = 0; i < hist_size ; i++){
        printf("%.4f ",hist1[i]);
        if(hist1[i]>max_color_v ){
            max_color_v  = hist1[i];
            hist1_max = i;
        }
    }
    printf("%d\n",hist1_max);
    //
    max_color_v = -1;
    for(int i = 0; i < hist_size ; i++){
        printf("%.4f ",hist2[i]);
        if(hist2[i]>max_color_v ){
            max_color_v  = hist2[i];
            hist2_max = i;
        }
    }
    printf("%d\n",hist2_max);    
    //printf("%d %d\n",max_color_idx1,max_color_idx2);    
    //
    max_color_v = -1;
    if (hist2_max == hist1_max){ // remove the background color
        for(int i=0;i<hist_size ; i++){ // find the second maximum
            if(i==hist1_max)
                continue;
            if(max_color_v<hist1[i]){
                max_color_v = hist1[i]; 
                hist1_max_new = i;            
            }
        }
    }else{
        hist1_max_new = hist1_max;
    }
    //printf("%s %d\n",sCTypes[hist1_max],hist1_max);
    //
    return hist1_max_new;
 }


int main(int argc, char** argv){    
	char *strFileImage;
    double cropPercentage=0.7;
    int num_color_blocks = 23;
    string buff;
    ifstream ifs("config");
    if(!ifs.fail()){
        getline(ifs,buff);
        num_color_blocks = atoi(buff.c_str());; // 0 black 1 white 2 gray
        cout << "Load config: num_color_block=" << num_color_blocks << endl;
    }else{
        cout << "No config: num_color_block=23" << endl;
    }
    int num_color_types = NUM_COLOR_TYPES;
    if (argc == 3){
		strFileImage = argv[1];	// Get image from first argument.
        cropPercentage = atof(argv[2]);
    }else{	
        cerr << "usage: MainColorDetection image_file center_percentage" << endl;
		cerr << "ERROR: No image was given on the command line!" << endl;
        for(int i=0;i<num_color_types ; i++){
            //printf("(%d)%s ",i,sCTypes[i]);            
            cerr << i << " " << sCTypes[i] << " ";
        }
        cerr << endl;
		return 1;
	}
	// Open the image, either as greyscale or color
	IplImage* imageIn = cvLoadImage(strFileImage, CV_LOAD_IMAGE_UNCHANGED);
	if (!imageIn) {
		cerr << "Couldn't load image file '" << strFileImage << "'" << endl;
		exit(1);
	}
	if (imageIn->nChannels != 3) {
		cerr << "Input image isn't a color RGB image!" << endl;
		exit(1);
	}
    //
    int hIn = imageIn->height;
    int wIn = imageIn->width;
    CvRect crop_rect;
    crop_rect.width = int(wIn * cropPercentage);
    crop_rect.height= int(hIn * cropPercentage);
    crop_rect.x = int((wIn - crop_rect.width)/2);
    crop_rect.y = int((hIn - crop_rect.height)/2);
    //
    IplImage *image_center = cvCreateImage(cvSize(crop_rect.width,crop_rect.height), 8, 3);
    cvSetImageROI(imageIn, crop_rect);
    cvCopy(imageIn, image_center, NULL);
    cvResetImageROI(imageIn);
    //
    double* colorHistogram1 = new double[num_color_types]; 
    double* colorHistogram2 = new double[num_color_types];
    for(int i = 0 ; i<num_color_types;i++){
        colorHistogram1[i]=0;
        colorHistogram2[i]=0;
    }
    double* colorHistogram_block1 = new double[num_color_blocks]; 
    double* colorHistogram_block2 = new double[num_color_blocks]; 
    for(int i = 0 ; i < num_color_blocks ; i++){
        colorHistogram_block1[i]=0;
        colorHistogram_block2[i]=0;
    }
    int max_color_idx1 = -1;
    int max_color_idx2 = -1;
    // color_type
    getColorHistogram(image_center , cvRect(0,0,0,0) , colorHistogram1 , num_color_types , 0); // the center image
    getColorHistogram(imageIn         , crop_rect         , colorHistogram2 , num_color_types , 0); // the border image
    max_color_idx1 = combineHistogram(num_color_types , colorHistogram1 , colorHistogram2);
    cout << sCTypes[max_color_idx1] << " " << max_color_idx1 << endl;
    // color_block
    getColorHistogram(image_center , cvRect(0,0,0,0) , colorHistogram_block1 , num_color_blocks , 1); // the center image
    getColorHistogram(imageIn         , crop_rect         , colorHistogram_block2 , num_color_blocks , 1); // the border image
    max_color_idx2 = combineHistogram(num_color_blocks , colorHistogram_block1 , colorHistogram_block2);
    cout <<  "Color Block Index " << max_color_idx2 << endl;
    //
    cvReleaseImage(&imageIn);
    cvReleaseImage(&image_center);
    delete []colorHistogram_block1;
    delete []colorHistogram_block2;
    delete []colorHistogram1;
    delete []colorHistogram2;
}

void getColorHistogram(IplImage* imageIn, const CvRect& croprect , double colorHistogram[] , int num_color_types , int type){
    int hIn = imageIn->height;
    int wIn = imageIn->width;

    // Create a HSV image showing the color types of the whole image, for debugging.
    IplImage *imageInHSV = cvCreateImage(cvGetSize(imageIn), 8, 3);
    cvCvtColor(imageIn, imageInHSV, CV_BGR2HSV);	// (note that OpenCV stores RGB images in B,G,R order.
    char *imOfsIn = imageInHSV->imageData;	// Pointer to the start of the input image HSV pixels.
    int rowSizeIn  = imageInHSV->widthStep;		// Size of row in bytes, including extra padding
    //    
    int ctype;
    for (int y=0; y<hIn; y++) {
        for (int x=0; x<wIn; x++) {
            // Get the HSV pixel components
            uchar H = *(uchar*)(imOfsIn + y*rowSizeIn + x*3 + 0);	// Hue
            uchar S = *(uchar*)(imOfsIn + y*rowSizeIn + x*3 + 1);	// Saturation
            uchar V = *(uchar*)(imOfsIn + y*rowSizeIn + x*3 + 2);	// Value (Brightness)
            // Determine what type of color the HSV pixel is.
            if(type==0){
#ifdef MY_DEBUG
                if(x==0&&y==0)
                    cout << "color_type" << endl;
#endif
                ctype = getPixelColorType(H, S, V);
            }else if(type==1){
#ifdef MY_DEBUG
                if(x==0&&y==0)
                    cout << "block_type" << endl;
#endif
                ctype = getPixelColorBlock(H , S , V , num_color_types);
            }else{
                cout<<"ERROR TYPE"<<endl;
                exit(0);
            }
            //Update Histogram
            if(croprect.height == 0 || croprect.width == 0){
                colorHistogram[ctype]++;
            }else{                
                if(!(x>croprect.x && x<croprect.x+croprect.width && y>croprect.y && y<croprect.y+croprect.height)){
                    colorHistogram[ctype]++;
                }
            }
        }
    }
    //    
    double max_color_v = -1;
    for(int i=0 ; i<num_color_types ; i++){
        colorHistogram[i] /= hIn * wIn;
        //if(max_color_v<colorHistogram[i]){
        //    max_color_v = colorHistogram[i]; 
        //    max_color_idx = i;            
        //}
        //printf("%s: %.4f\n",sCTypes[i],colorHistogram[i]);
    }
    //printf("Main Color: %s\n",sCTypes[max_color_idx]);
    //
	cvReleaseImage(&imageInHSV);    
}