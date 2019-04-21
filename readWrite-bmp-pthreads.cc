
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <time.h>
#include <omp.h>
#include <chrono>

#define WIDTH	3
#define HEIGHT	3
#define DEBUG	1
#define KERNEL_SIZE	3

int threshold = 20;

using namespace std;
using namespace std::chrono;

#pragma pack(1)
typedef struct {
	char id[2];
	int file_size;
	int reserved;
	int offset;

} header_type;

#pragma pack(1)
typedef struct {
	int header_size;
	int width;
	int height;
	unsigned short int color_planes;
	unsigned short int color_depth;
	unsigned int compression;
	int image_size;
	int xresolution;
	int yresolution;
	int num_colors;
	int num_important_colors;
} information_type;

int globalPrintCount = 0; //this will keep count of the times printArray is called

void printImageArray(vector<vector<int> > &image) 
{
	int imageRows = image.size();
	int imageColsIndex = image[0].size();
	cout << "print time: " << globalPrintCount << endl;
	globalPrintCount++;

	for (int i = 0; i < imageRows; i++) {
		cout << "Row " << i << "[ ";
		for (int j = 0; j < imageColsIndex; j++) {
			cout << image[i][j] << "\t";
		}
		cout << " ]" << endl;
	}
	cout << endl << endl;
}

unsigned char tempData[3];
int row, col, row_bytes, padding;
vector <vector <int> > data, newData;

ifstream imageFile;
ofstream newImageFile;

vector< vector<int> > prepareKernel(int x, int y/*, vector< vector<int> > &image*/)
{

	int imageRowsIndex = data.size() - 1;
	int imageColsIndex = data[0].size() - 1;

	vector<vector<int> > kernel(KERNEL_SIZE, vector<int>(KERNEL_SIZE, 0));

	if (x > 0 && y > 0 && x < imageRowsIndex && y < imageColsIndex) {
		for (int i = 0; i < KERNEL_SIZE; i++) {
			for (int j = 0; j < KERNEL_SIZE; j++) {
				kernel[i][j] = data[(x - 1) + i][(y - 1) + j];
			}

		}
	} else if (x == 0 && y > 0 && x < imageRowsIndex && y < imageColsIndex) {
		for (int i = 1; i < KERNEL_SIZE; i++) {
			for (int j = 0; j < KERNEL_SIZE; j++) {
				kernel[i][j] = data[(x - 1) + i][(y - 1) + j];
			}
		}
	} else if (x == 0 && y == 0 && x < imageRowsIndex && y < imageColsIndex) {
		for (int i = 1; i < KERNEL_SIZE; i++) {
			for (int j = 1; j < KERNEL_SIZE; j++) {
				kernel[i][j] = data[(x - 1) + i][(y - 1) + j];
			}
		}
	} else if (x > 0 && y == 0 && x < imageRowsIndex && y < imageColsIndex) {
		for (int i = 0; i < KERNEL_SIZE; i++) {
			for (int j = 1; j < KERNEL_SIZE; j++) {
				kernel[i][j] = data[(x - 1) + i][(y - 1) + j];
			}
		}
	} /*else if (x == imageRowsIndex && y < imageColsIndex) {
		for (int i = 0; i < KERNEL_SIZE - 1; i++) {
			for (int j = 0; j < KERNEL_SIZE; j++) {
				kernel[i][j] = data[(x - 1) + i][(y - 1) + j];
			}
		}
	}*/ else if (x == imageRowsIndex && y == 0) {
		for (int i = 0; i < KERNEL_SIZE - 1; i++) {
			for (int j = 1; j < KERNEL_SIZE; j++) {
				kernel[i][j] = data[(x - 1) + i][(y - 1) + j];
			}
		}
	} else if (x == imageRowsIndex && y == imageColsIndex) {
		for (int i = 0; i < KERNEL_SIZE - 1; i++) {
			for (int j = 0; j < KERNEL_SIZE - 1; j++) {
				kernel[i][j] = data[(x - 1) + i][(y - 1) + j];
			}
		}
	} else if (x == 0 && y == imageColsIndex) {
		for (int i = 1; i < KERNEL_SIZE; i++) {
			for (int j = 0; j < KERNEL_SIZE - 1; j++) {
				kernel[i][j] = data[(x - 1) + i][(y - 1) + j];
			}
		}
	} else if (x < imageRowsIndex && y == imageColsIndex) {
		for (int i = 0; i < KERNEL_SIZE; i++) {
			for (int j = 0; j < KERNEL_SIZE - 1; j++) {
				kernel[i][j] = data[(x - 1) + i][(y - 1) + j];
			}
		}
	}

	return kernel;
}



int sobelFilter (int x, int y)
{
	//prepare data;
	vector< vector<int> > dataToFilter = prepareKernel(x,y);

	// printImageArray(dataToFilter);

	//myfilter's algorithm is quite simple, we just add
	//the derivatives in x and y (called gradient magnitude)
	int x0,x1,x2,x3,x5,x6,x7,x8;
	x0 = dataToFilter[0][0];
	x1 = dataToFilter[0][1];
	x2 = dataToFilter[0][2];
	x3 = dataToFilter[1][0];
	x5 = dataToFilter[1][2];
	x6 = dataToFilter[2][0];
	x7 = dataToFilter[2][1];
	x8 = dataToFilter[2][2];

	int dfdy = (x0 + 2*x1 +x2) - (x6 + 2*x7 + x8);
	int dfdx = (x2 + 2*x5 +x8) - (x0 + 2*x3 + x6);

	int gradient = abs(dfdy) + abs(dfdx);

	//the actual filter:
	if(gradient < threshold)
	{
		return 0;
	} else {
		return 255;
	}
}

struct index_limits {
	int lower_index;
	int higher_index;
	int threadId;
};

void *ProcessData(void *arg) {
	clock_t t = clock();
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	double startTime = omp_get_wtime();

	struct index_limits *limits = (struct index_limits *)arg;

	int lower_index = limits->lower_index;
	int higher_index = limits->higher_index;

	if((size_t)lower_index > data.size()){
		lower_index = data.size();
	}

	if((size_t)higher_index > data.size()){
		higher_index = data.size();
	}
	printf("Thread id %d lower_index:%d higher_index%d\n",limits->threadId, lower_index, higher_index);
	for (int row=lower_index; row < higher_index; row++) {
		for (int col=0; (size_t)col < data[0].size(); col++) {
			newData[row][col] = sobelFilter(row,col);
		}
	}

	t = clock() - t;
	double stopTime = omp_get_wtime();
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	float secsElapsed_clock = (float)t/CLOCKS_PER_SEC;
	double secsElapsed_omp = stopTime - startTime;
	auto secsElapsed_high = duration_cast<microseconds>(t2-t1).count();

	printf("Thread ID: %d - Using clock Elapsed time: %f \n", limits->threadId, secsElapsed_clock);
	printf("Thread ID: %d - Using omp_get_wtime Elapsed time: %f \n", limits->threadId, secsElapsed_omp);
	printf("Thread ID: %d - Using HighRes Elapsed time: %ld \n", limits->threadId, secsElapsed_high);

	return NULL;
}

int main(int argc, char* argv[])
{
	int numThreads = 1;
	string imageFileName, newImageFileName;
	header_type header;
	information_type information;

	string strThreshold, strNumThreads;
	
	if ( argc == 1 ) {
		// prepare files
		cout << "Original imagefile? ";
		cin >> imageFileName;
		
		cout << "New imagefile name? ";
		cin >> newImageFileName;

		cout << "Threshold? ";
		cin >> strThreshold;

		//Request number of threads
		cout << "Number of threads? ";
		cin >> strNumThreads;
	} else if ( argc==5 ) {
		imageFileName = argv[1];
		newImageFileName = argv[2];
		strThreshold=argv[3];
		strNumThreads=argv[4];
	} else {
		cout << "Usage: "<< endl<<"./readWrite-bmp-threaded [<original image path> <new image name> <threshold> <# of threads>]"<<endl;
		return 1;
	}
		

	imageFile.open(imageFileName.c_str(), ios::binary);

        if (!imageFile) {
                cerr << "file not found" << endl;
                exit(-1);
        }
	
	newImageFile.open(newImageFileName.c_str(), ios::binary);

	// read file header
	imageFile.read((char *) &header, sizeof(header_type));
	if (header.id[0] != 'B' || header.id[1] != 'M') {
		cerr << "Does not appear to be a .bmp file.  Goodbye." << endl;
		exit(-1);
	}

	threshold = atoi(strThreshold.c_str());

	numThreads = atoi(strNumThreads.c_str());

	// read/compute image information
	imageFile.read ((char *) &information, sizeof(information_type));
	row_bytes = information.width * 3;

	padding = row_bytes % 4;
	if (padding)
		padding = 4 - padding;


	// extract image data, initialize vectors
	for (row=0; row < information.height; row++) {
		data.push_back (vector <int>());
		for (col=0; col < information.width; col++) {
			imageFile.read ((char *) tempData, 3 * sizeof(unsigned char));
			data[row].push_back ((int) tempData[0]);
		}
		if (padding)
			imageFile.read ((char *) tempData, padding * sizeof(unsigned char));
	}
	cout << imageFileName << ": " << information.width << " x " << information.height << endl;

	for(int i = 0; i < information.height; i++){
		newData.push_back(vector<int>(information.width, 0));
	}

	clock_t clock_time = clock();
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	double startTime = omp_get_wtime();

	pthread_t threads[numThreads];
	int rc;
	long t;

	struct index_limits args[numThreads] ;

	float cs_float = floor(information.height/numThreads);
	int chunck_size = (cs_float >= 0) ? (int)(cs_float + 0.5) : (int)(cs_float - 0.5);

	for(t=0; t < numThreads; t++){

		printf("In main: creating thread %ld\n", t);
		args[t].lower_index = chunck_size*t;

		if(t != numThreads -1){
			args[t].higher_index = args[t].lower_index + chunck_size;
		} else {
			args[t].higher_index = information.height;
		}
		args[t].threadId = t;
		rc = pthread_create(&threads[t], NULL, ProcessData, (void *)&args[t]);
		if (rc){
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}

	void* status;
	for(int t = 0; t < numThreads; t++)
	{
		rc = pthread_join(threads[t],&status);
		if(rc){
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
	}



	clock_time = clock() - clock_time;
	double stopTime = omp_get_wtime();
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	float secsElapsed_clock = (float)clock_time/CLOCKS_PER_SEC;
	double secsElapsed_omp = stopTime - startTime;
	auto secsElapsed_high = duration_cast<microseconds>(t2-t1).count();

	//printf("General threaded processing time - Using clock Elapsed time: %f \n", (secsElapsed_clock));
	printf("General threaded processing time - Using omp_get_wtime Elapsed time: %f \n", secsElapsed_omp);
	//printf("General threaded processing time - Using HighRes Elapsed time: %ld \n", secsElapsed_high);

	// write header to new image file
	newImageFile.write ((char *) &header, sizeof(header_type));
	newImageFile.write ((char *) &information, sizeof(information_type));

	// write new image data to new image file
	for (row=0; row < information.height; row++) {
		for (col=0; col < information.width; col++) {
			tempData[0] = (unsigned char) newData[row][col];
			tempData[1] = (unsigned char) newData[row][col];
			tempData[2] = (unsigned char) newData[row][col];
			newImageFile.write ((char *) tempData, 3 * sizeof(unsigned char));
		}

		if (padding) {
			tempData[0] = 0;
			tempData[1] = 0;
			tempData[2] = 0;
			newImageFile.write ((char *) tempData, padding * sizeof(unsigned char));
		}
	}
	cout << newImageFileName << " done." << endl;
	imageFile.close();
	newImageFile.close();

	/* Last thing that main() should do */
	pthread_exit(NULL);
	return 0;
}

