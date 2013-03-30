//uploaded to google code with message
//#define SHOW_DEBUG_IMAGE
#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
//#include <stdio.h>
//#include <tchar.h>

#include <cstdio>	// Used for "printf"
#include <string>	// Used for C++ strings
#include <iostream>	// Used for C++ cout print statements
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
//#include <cmath>	// Used to calculate square-root for statistics

// Include OpenCV libraries
#include <cv.h>
#include <cvaux.h>
#include <cxcore.h>
#include <highgui.h>
#include <ml.h>
#include <math.h>
using namespace std;

// Various color types for detected shirt colors.
enum                             {cBLACK=0,cWHITE, cGREY, cRED, cORANGE, cBrown, cYELLOW, cGREEN, cCyan, cBLUE, cPURPLE, cPINK,  NUM_COLOR_TYPES};
char* sCTypes[NUM_COLOR_TYPES] = {"Black"   , "White"   , "Grey"    , "Red" , "Orange"  , "Brown"   , "Yellow"  , "Green"   , "Cyan"    , "Blue" , "Purple" , "Pink"};
uchar cCTHue[NUM_COLOR_TYPES]  = {0         , 0         , 0         , 0     , 20        , 20        , 30        , 55        , 85        , 115    , 138      , 161};
uchar cCTSat[NUM_COLOR_TYPES]  = {0         , 0         , 0         , 255   , 255       , 255       , 255       , 255       , 255       , 255    , 255      , 255};
uchar cCTVal[NUM_COLOR_TYPES]  = {0         , 255       , 120       , 255   , 255       , 165       , 255       , 255       , 255       , 255    , 255      , 255};

typedef struct ColorHistogramItem{
    int colorIndex;
    double histv;    
}ColorHistogramItem;

CvKNearest *knn = NULL;
CvMat *colorCenterData = NULL;
CvMat *colorCenterClass = NULL;

int getPixelColorType(int H, int S, int V);
// for S, three range 1/4 - 2/4 - 3/4 - 4/4
// for V, four range 0/4 - 1/4 - 2/4 - 3/4 - 4/4
// for each H range, S-V id is : 1 2 3 / 4 5 6 / 7 8 9 / 10 11 12
int getPixelColorBlock(int H , int S , int V , int num_color_types);
int getPixelColorCluster(int H, int S, int V, int R , int G , int B);
int getPixelColorDistance(int H, int S, int V, int R , int G , int B, CvMat* colorCenterData);
void getColorHistogram(IplImage* img, const CvRect& croprect, ColorHistogramItem colorHistogram[], int num_color_types, int type=0);
void getColorHistogramDistance(IplImage* img, const CvRect& croprect, ColorHistogramItem colorHistogram[], int num_color_types, CvMat* colorCenterData);
// hist1 -> center hist2 -> background
//void outputMainColor(ColorHistogramItem hist_center[], ColorHistogramItem hist_board[], int hist_size);
void outputMainColor(ColorHistogramItem* hist_center, ColorHistogramItem* hist_board, int hist_size);

int main(int argc, char** argv){    
	char *strFileImage;
    double cropPercentage=0.7;
    int num_color_blocks = 4 * 16 + 3; // default H has 4 blocks
    string buff;
    //
    ifstream ifs("config");
    if(!ifs.fail()){
        getline(ifs,buff);
        num_color_blocks = atoi(buff.c_str()) * 16 + 3; // 0 black 1 white 2 gray
        cout << "Load config: num_color_block=" << num_color_blocks << "(0 black 1 white 2 gray + H_blocks * S_blocks(4) * V_blocks(4) )"<<endl;
    }else{
        cout << "No config: num_color_block = 4 * 16 + 3 = 67" << endl;
    }
    ifs.close();
    //
    double rfv,gfv,bfv;
    int cluster_NUM, cluster_count=-1;
    ifs.open("/home/dayong/WORK/maincolor/colorlist_rgb.txt", ifstream::in);    
    if(!ifs.fail()){
        getline(ifs,buff);
        cluster_NUM = atoi(buff.c_str());
        colorCenterData = cvCreateMat( cluster_NUM , 3 , CV_32FC1 );
        colorCenterClass = cvCreateMat( cluster_NUM , 1 , CV_32FC1 );    
        //
        getline(ifs,buff);
        do{            
            cluster_count++;
            sscanf(buff.c_str(),"%lf %lf %lf",&rfv,&gfv,&bfv);            
            cvmSet(colorCenterData, cluster_count, 0, rfv);
            cvmSet(colorCenterData, cluster_count, 1, gfv);
            cvmSet(colorCenterData, cluster_count, 2, bfv);
            cvmSet(colorCenterClass, cluster_count, 0, cluster_count);
            //cout << cluster_count << " " << rfv <<" " << gfv << " " << bfv << endl;
            getline(ifs,buff);
        }while(!ifs.eof());
        //knn = new CvKNearest( colorCenterData, colorCenterClass, 0, false, 1 );
        //CvMat* nearests = cvCreateMat(1 , 1 , CV_32FC1);
    }else{
        cout << "Fail to open colorlist_rgb.txt" << endl;
        exit(0);
    }
    ifs.close();
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
        if (imageIn->nChannels != 1) {
    		cerr << "Input image isn't a color RGB image or black-white image: " << imageIn->nChannels << endl;
    		exit(1);
        }else{
            IplImage* imageIn_tmp = imageIn;
            imageIn = cvCreateImage(cvGetSize(imageIn_tmp),imageIn_tmp->depth,3);
            cvMerge(imageIn_tmp,NULL,NULL,NULL,imageIn);
            cvReleaseImage(&imageIn_tmp);
        }   
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
    
    // color_type
    ColorHistogramItem* colorHistogram1 = new ColorHistogramItem[num_color_types]; 
    ColorHistogramItem* colorHistogram2 = new ColorHistogramItem[num_color_types];
    for(int i = 0 ; i<num_color_types;i++){        
        colorHistogram1[i].colorIndex = i;
        colorHistogram2[i].colorIndex = i;
        colorHistogram1[i].histv=0;
        colorHistogram2[i].histv=0;
    }    
    getColorHistogram(image_center , cvRect(0,0,0,0) , colorHistogram1 , num_color_types , 0); // the center image
    getColorHistogram(imageIn         , crop_rect         , colorHistogram2 , num_color_types , 0); // the border image
    outputMainColor(colorHistogram1, colorHistogram2, num_color_types);    
    
    // color block
    //ColorHistogramItem* colorHistogram_block1 = new ColorHistogramItem[num_color_blocks]; 
    //ColorHistogramItem* colorHistogram_block2 = new ColorHistogramItem[num_color_blocks]; 
    //for(int i = 0 ; i < num_color_blocks ; i++){
    //    colorHistogram_block1[i].colorIndex = i;
    //    colorHistogram_block2[i].colorIndex = i;
    //    colorHistogram_block1[i].histv=0;
    //    colorHistogram_block2[i].histv=0;
    //}
    //int max_color_idx1, max_color_idx2;    
    // color_block
    //getColorHistogram(image_center , cvRect(0,0,0,0) , colorHistogram_block1 , num_color_blocks , 1); // the center image
    //getColorHistogram(imageIn         , crop_rect        , colorHistogram_block2 , num_color_blocks , 1); // the border image
    //max_color_idx2 = combineHistogram(num_color_blocks , colorHistogram_block1 , colorHistogram_block2);
    //cout <<  "Color Block Index " << max_color_idx2 << endl;
    
    // color_type
    ColorHistogramItem* colorHistogram_cluster1 = new ColorHistogramItem[cluster_NUM+3]; 
    ColorHistogramItem* colorHistogram_cluster2 = new ColorHistogramItem[cluster_NUM+3];
    for(int i = 0 ; i<cluster_NUM+3 ; i++){        
        colorHistogram_cluster1[i].colorIndex = i;
        colorHistogram_cluster2[i].colorIndex = i;
        colorHistogram_cluster1[i].histv=0;
        colorHistogram_cluster2[i].histv=0;
    }    
    getColorHistogramDistance(image_center , cvRect(0,0,0,0) , colorHistogram_cluster1 , cluster_NUM+3 , colorCenterData); // the center image
    getColorHistogramDistance(imageIn         , crop_rect         , colorHistogram_cluster2 , cluster_NUM+3 , colorCenterData); // the border image
    outputMainColor(colorHistogram_cluster1, colorHistogram_cluster2, cluster_NUM+3);    
    
    // release all the sources
    cvReleaseImage(&imageIn);
    cvReleaseImage(&image_center);
    delete []colorHistogram1;
    delete []colorHistogram2;
    //delete []colorHistogram_block1;
    //delete []colorHistogram_block2;
    delete []colorHistogram_cluster1;
    delete []colorHistogram_cluster2;    
    delete knn;
    cvReleaseMat(&colorCenterClass);
    cvReleaseMat(&colorCenterData);
}

static int colorHistogramItemCompare(const void* a, const void* b){
    ColorHistogramItem* item1 = (ColorHistogramItem*)a;
    ColorHistogramItem* item2 = (ColorHistogramItem*)b;
    return item1->histv  < item2->histv ? 1 : -1;
}

int getPixelColorType(int H, int S, int V)
{
	int color;
	if (V < 75) // check it is black
		color = cBLACK;
	else if (V > 175 && S < 47) // original 47
		color = cWHITE;
	else if (V < 175 && S < 63) // original 63
		color = cGREY;
	else {	// Is a color
		if (H < 14)
			color = cRED;
        else if (H < 25){
            if(V<165){
			    color = cBrown;
            }else{
                color = cORANGE;
            }
        }else if (H < 34)
			color = cYELLOW;
		else if (H < 73)
			color = cGREEN;
		else if (H < 102)
			color = cCyan;
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

int getPixelColorBlock(int H , int S , int V , int num_color_types){
	int color, h_index, s_index, v_index, sv_table_index;    
    int color_block = num_color_types  - 3; // 0 black 1 white 2 gray
    int h_blocks  = color_block / 16;
    int h_block_width = 180 / h_blocks;    
	if (V < 75) // if v < 1/4 * 255, then it is black
		color = 0; // black
	else if (V > 180 && S < 47) // if v >  (5/8 -> 6/8) * 255 and s < 1/8 * 255, then it is white
		color = 1; // white
	else if (V < 180 && S < 63) // if v <  (5/8 -> 6/8) * 255 and s < 1/4 * 255, then it is white
		color = 2; // gray
	else {	// Is a color
        h_index = -1;
        int H_ori = H;
        while(H>=0){
            h_index++;
            H -= h_block_width;
        }
        if(h_index >= h_blocks){ // in case H = 180
            h_index = h_blocks-1;
        }
        v_index = floor( V / 64.0 );
        if(v_index<0 || v_index>4){
            cerr << "ERROR: v_index < 0 || v_index > 4" << endl;
            exit(0);
        }
        s_index = floor( S / 64.0 );
        if(s_index < 0 || s_index > 4){
            cerr << "ERROR: s_index < 0 || s_index > 4" << endl;            
            exit(0);
        }
        sv_table_index = v_index * 4 + s_index;
        color = h_index * 16 + sv_table_index + 3;
        if(color<0 || color >= num_color_types){
            cerr << "ERROR: color<0 || color >=" << num_color_types << " " << H_ori << " " << h_index << " " << s_index << " " << v_index << endl;            
            exit(0);
        }
        //color +=3;
    }
    return color;
}

int getPixelColorCluster(int H , int S , int V , int R , int G , int B){
    int K = 1, color;
    float _sample[3], response;    
    CvMat sample   =  cvMat( 1, 3, CV_32FC1, _sample );        
	if (V < 75) // check it is black
		color = 0;
	else if (V > 175 && S < 47) // original 47
		color = 1;
	else if (V < 175 && S < 63) // original 63
		color = 2;
	else {	// Is a color
        sample.data.fl[0] = (float)R/255.0;
        sample.data.fl[1] = (float)G/255.0;
        sample.data.fl[2] = (float)B/255.0;
        response = knn->find_nearest(&sample,K,0,0,NULL,0);
        color = response + 3;
    }
    return color;
}

static inline double getDistance(double R, double G, double B, double cR, double cG, double cB){
    long mR = (R + cR)/2;
    double mR1 = mR >> 8;
    double mR2 = (1+mR) >> 8;
    double dR = R - cR;
    double dG = G - cG;
    double dB = B - cB;
    //double dis = sqrt( ((512+mR)*dR*dR)>>8 + 4 * dG * dG + ((767-mR)*dB*dB)>>8 );
    double dis = ((2+mR1)*dR*dR) + 4 * dG * dG + ((3-mR2)*dB*dB);
    return dis;
}

int getPixelColorDistance(int H, int S, int V, int R , int G , int B, CvMat* colorCenterData){
    int K = 1, color;
    float _sample[3], response;
    double cR,cG,cB,dis,min_dis;
    int min_dis_idx;
    CvMat sample   =  cvMat( 1, 3, CV_32FC1, _sample );        
	if (V < 75) // check it is black
		color = 0;
	else if (V > 175 && S < 47) // original 47
		color = 1;
	else if (V < 175 && S < 63) // original 63
		color = 2;
	else {	// Is a color
        for(int i = 0 ; i<colorCenterData->rows ; i++){
            cR = floor(cvGetReal2D(colorCenterData,i,0)*255.0);
            cG = floor(cvGetReal2D(colorCenterData,i,1)*255.0);
            cB = floor(cvGetReal2D(colorCenterData,i,2)*255.0);
            dis = getDistance((double)R, (double)G, (double)B, cR, cG, cB);
            if (i==0){
                min_dis = dis;
                min_dis_idx = 0;
            }else{
                if(dis<min_dis){
                    min_dis = dis;
                    min_dis_idx = i;
                }
            }
        }
        color = min_dis_idx + 3;
    }
    return color;
}

void getColorHistogram(IplImage* imageIn, const CvRect& croprect , ColorHistogramItem colorHistogram[] , int num_color_types , int type){
    int hIn = imageIn->height;
    int wIn = imageIn->width;

    // Create a HSV image showing the color types of the whole image, for debugging.
    IplImage *imageInHSV = cvCreateImage(cvGetSize(imageIn), 8, 3);
    cvCvtColor(imageIn, imageInHSV, CV_BGR2HSV);	// (note that OpenCV stores RGB images in B,G,R order.
    char *imOfsIn1 = imageIn->imageData;
    int rowSizeIn1 = imageIn->widthStep;
    char *imOfsIn2 = imageInHSV->imageData;	// Pointer to the start of the input image HSV pixels.
    int rowSizeIn2  = imageInHSV->widthStep;		// Size of row in bytes, including extra padding
#ifdef _DEBUG
    IplImage *imageResult = cvCreateImage(cvGetSize(imageIn), 8, 3);
#endif
    //    
    int ctype;
    for (int y=0; y<hIn; y++) {
        for (int x=0; x<wIn; x++) {
            // Get the RGB pixel components
            uchar B  = *(uchar*)(imOfsIn1 + y*rowSizeIn1 + x*3 + 0);	// B
            uchar G  = *(uchar*)(imOfsIn1 + y*rowSizeIn1 + x*3 + 1);	// G
            uchar R  = *(uchar*)(imOfsIn1 + y*rowSizeIn1 + x*3 + 2);	// R
            // Get the HSV pixel components
            uchar H = *(uchar*)(imOfsIn2 + y*rowSizeIn2 + x*3 + 0);	// Hue
            uchar S  = *(uchar*)(imOfsIn2 + y*rowSizeIn2 + x*3 + 1);	// Saturation
            uchar V = *(uchar*)(imOfsIn2 + y*rowSizeIn2 + x*3 + 2);	// Value (Brightness)
            // Determine what type of color the HSV pixel is.
            if(type==0){
                ctype = getPixelColorType(H, S, V);
            }else if(type==1){
                ctype = getPixelColorBlock(H , S , V , num_color_types);
            } else if(type==2){
                ctype = getPixelColorCluster(H, S, V, R, G, B);
#ifdef _DEBUG
                if(ctype == 0){
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 2) = 0;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 1) = 0;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 0) = 0;
                }else if(ctype == 1){
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 2) = 255;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 1) = 255;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 0) = 255;
                }else if(ctype == 2){
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 2) = 123;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 1) = 123;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 0) = 123;
                }else if(ctype>=3){
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 2) = floor(cvGetReal2D(colorCenterData,ctype-3 , 0)*255);
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 1) = floor(cvGetReal2D(colorCenterData,ctype-3 , 1)*255);
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 0) = floor(cvGetReal2D(colorCenterData,ctype-3 , 2)*255);
                }
#endif
            } else {
                cout<<"ERROR TYPE"<<endl;
                exit(0);
            }
            //Update Histogram
            if(ctype<0 || ctype>=num_color_types){
                cout << "ERROR: (ctype<0 || ctype>=num_color_types)" << ctype << endl;
            }
            if(croprect.height == 0 || croprect.width == 0){
                colorHistogram[ctype].histv= colorHistogram[ctype].histv + 1.0;
            }else{                
                if(!(x>croprect.x && x<croprect.x+croprect.width && y>croprect.y && y<croprect.y+croprect.height)){
                    colorHistogram[ctype].histv = colorHistogram[ctype].histv + 1.0;
                }
            }
        }
    }
    //
#ifdef _DEBUG
    cvNamedWindow("DEBUG", CV_WINDOW_AUTOSIZE); 
    cvShowImage("DEBUG", imageResult);
    cvWaitKey(0);
#endif
    //
    qsort(colorHistogram,num_color_types,sizeof(ColorHistogramItem),colorHistogramItemCompare);
    //    
    double max_color_v = -1;
    double hist_sum = 0;    
    for(int i=0 ; i<num_color_types ; i++){
        hist_sum += colorHistogram[i].histv;
    }    
    for(int i=0 ; i<num_color_types ; i++){
        colorHistogram[i].histv = colorHistogram[i].histv /  hist_sum;
    }
#ifdef _DEBUG
    for(int i=0 ; i<num_color_types ; i++){
        printf("(%d)%.4lf ",colorHistogram[i].colorIndex,colorHistogram[i].histv);
    }
    printf("\n");   
#endif
    cvReleaseImage(&imageInHSV);    
#ifdef _DEBUG
    cvReleaseImage(&imageResult);
#endif
}

int old_combineHistogram(int hist_size, double hist1[], double hist2[]){
    int hist1_max, hist2_max, hist1_max_new;
    double max_color_v  = -1;
   //
    max_color_v  = -1;
    for(int i = 0; i < hist_size ; i++){
        printf("%.4lf ",hist1[i]);
        if(hist1[i]>max_color_v ){
            max_color_v  = hist1[i];
            hist1_max = i;
        }
    }
    printf("%d\n",hist1_max);
    //
    max_color_v = -1;
    for(int i = 0; i < hist_size ; i++){
        printf("%.4lf ",hist2[i]);
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

void outputMainColor(ColorHistogramItem hist_center[], ColorHistogramItem hist_board[], int hist_size){    
    int i,topN = 2;
    double sumv = 0.0;
    if(hist_center[0].colorIndex == hist_board[0].colorIndex && hist_center[0].histv <= 0.98){ // the first color is background color and it is the only color in the center        
        for(i = 1 ; i < hist_size ; i++){
            sumv += hist_center[i].histv;
        }
        for(i = 1 ; i < 1 + topN ; i++){
            printf("%d %.4lf ",hist_center[i].colorIndex, hist_center[i].histv / sumv);
        }
        printf("\n");
    }else{
        if(hist_center[1].colorIndex == hist_board[0].colorIndex && hist_center[1].histv >= 0.3){ // the second color is background color and it has a large portion
            for(i = 0 ; i < hist_size ; i++){
                if(i==1)
                    continue;
                sumv += hist_center[i].histv;
            }            
            for(i = 0 ; i < 1 + topN ; i++){
                if(i==1)
                    continue;
                printf("%d %.4lf ",hist_center[i].colorIndex, hist_center[i].histv / sumv);
            }
            printf("\n");        
        }else{
             for(i = 0 ; i < topN ; i++){
                printf("%d %.4lf ",hist_center[i].colorIndex, hist_center[i].histv);
            }
            printf("\n");
        }
    }
 }


void getColorHistogramDistance(IplImage* imageIn, const CvRect& croprect , ColorHistogramItem colorHistogram[] , int num_color_types , CvMat* colorCenterData){
    int hIn = imageIn->height;
    int wIn = imageIn->width;

    // Create a HSV image showing the color types of the whole image, for debugging.
    IplImage *imageInHSV = cvCreateImage(cvGetSize(imageIn), 8, 3);
    cvCvtColor(imageIn, imageInHSV, CV_BGR2HSV);	// (note that OpenCV stores RGB images in B,G,R order.
    char *imOfsIn1 = imageIn->imageData;
    int rowSizeIn1 = imageIn->widthStep;
    char *imOfsIn2 = imageInHSV->imageData;	// Pointer to the start of the input image HSV pixels.
    int rowSizeIn2  = imageInHSV->widthStep;		// Size of row in bytes, including extra padding
#ifdef _DEBUG
    IplImage *imageResult = cvCreateImage(cvGetSize(imageIn), 8, 3);
#endif
    //    
    int ctype;
    for (int y=0; y<hIn; y++) {
        for (int x=0; x<wIn; x++) {
            // Get the RGB pixel components
            uchar B  = *(uchar*)(imOfsIn1 + y*rowSizeIn1 + x*3 + 0);	// B
            uchar G  = *(uchar*)(imOfsIn1 + y*rowSizeIn1 + x*3 + 1);	// G
            uchar R  = *(uchar*)(imOfsIn1 + y*rowSizeIn1 + x*3 + 2);	// R
            // Get the HSV pixel components
            uchar H = *(uchar*)(imOfsIn2 + y*rowSizeIn2 + x*3 + 0);	// Hue
            uchar S  = *(uchar*)(imOfsIn2 + y*rowSizeIn2 + x*3 + 1);	// Saturation
            uchar V = *(uchar*)(imOfsIn2 + y*rowSizeIn2 + x*3 + 2);	// Value (Brightness)
            // Determine what type of color the HSV pixel is.
               ctype = getPixelColorDistance(H, S, V, R, G, B, colorCenterData);
#ifdef _DEBUG
                if(ctype == 0){
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 2) = 0;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 1) = 0;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 0) = 0;
                }else if(ctype == 1){
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 2) = 255;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 1) = 255;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 0) = 255;
                }else if(ctype == 2){
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 2) = 123;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 1) = 123;
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 0) = 123;
                }else if(ctype>=3){
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 2) = floor(cvGetReal2D(colorCenterData,ctype-3 , 0)*255);
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 1) = floor(cvGetReal2D(colorCenterData,ctype-3 , 1)*255);
                    *(uchar*)(imageResult->imageData + y*rowSizeIn1 + x*3 + 0) = floor(cvGetReal2D(colorCenterData,ctype-3 , 2)*255);
                }         
#endif
            //Update Histogram
            if(ctype<0 || ctype>=num_color_types){
                cout << "ERROR: (ctype<0 || ctype>=num_color_types)" << ctype << endl;
            }
            if(croprect.height == 0 || croprect.width == 0){
                colorHistogram[ctype].histv= colorHistogram[ctype].histv + 1.0;
            }else{                
                if(!(x>croprect.x && x<croprect.x+croprect.width && y>croprect.y && y<croprect.y+croprect.height)){
                    colorHistogram[ctype].histv = colorHistogram[ctype].histv + 1.0;
                }
            }
        }
    }
    //
#ifdef _DEBUG
    cvNamedWindow("DEBUG", CV_WINDOW_AUTOSIZE); 
    cvShowImage("DEBUG", imageResult);
    cvWaitKey(0);
#endif
    //
    qsort(colorHistogram,num_color_types,sizeof(ColorHistogramItem),colorHistogramItemCompare);
    //    
    double max_color_v = -1;
    double hist_sum = 0;    
    for(int i=0 ; i<num_color_types ; i++){
        hist_sum += colorHistogram[i].histv;
    }    
    for(int i=0 ; i<num_color_types ; i++){
        colorHistogram[i].histv = colorHistogram[i].histv /  hist_sum;
    }
#ifdef _DEBUG
    for(int i=0 ; i<num_color_types ; i++){
        printf("(%d)%.4lf ",colorHistogram[i].colorIndex,colorHistogram[i].histv);
    }
    printf("\n");   
#endif
    cvReleaseImage(&imageInHSV);    
#ifdef _DEBUG
    cvReleaseImage(&imageResult);
#endif
}
