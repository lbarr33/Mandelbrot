#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
////////////////////////////////////////////////////////////////////////////////

const int MAX_WIDTH_HEIGHT = 30000;
const int HUE_PER_ITERATION = 5;
const bool DRAW_ON_KEY = true;
const int CHUNK_SIZE = 1;
const int PROCNUM = 32;
////////////////////////////////////////////////////////////////////////////////
double When();
int rank;
MPI_Status status;

class State {
    public:
        double centerX;
        double centerY;
        double zoom;
        int maxIterations;
        int w;
        int h;
        State() {
            //centerX = -.75;
            //centerY = 0;
            centerX = -1.186340599860225;
            centerY = -0.303652988644423;
            zoom = 1;
            maxIterations = 100;
            w = 28000;
            h = 28000;
        }
};

////////////////////////////////////////////////////////////////////////////////

float iterationsToEscape(double x, double y, int maxIterations) {
    double tempa;
    double a = 0;
    double b = 0;
    for (int i = 0 ; i < maxIterations ; i++) {
        tempa = a*a - b*b + x;
        b = 2*a*b + y;
        a = tempa;
        if (a*a+b*b > 64) {
            // return i; // discrete
            return i - log(sqrt(a*a+b*b))/log(8); //continuous
        }
    }
    return -1;
}

int hue2rgb(float t){
    while (t>360) {
        t -= 360;
    }
    if (t < 60) return 255.*t/60.;
    if (t < 180) return 255;
    if (t < 240) return 255. * (4. - t/60.);
    return 0;
}

void writeImage(unsigned char *img, int w, int h) {
    long long filesize = 54 + 3*(long long)w*(long long)h;
    unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
 unsigned char bmppad[3] = {0,0,0};

    bmpfileheader[ 2] = (unsigned char)(filesize    );
    bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
    bmpfileheader[ 4] = (unsigned char)(filesize>>16);
    bmpfileheader[ 5] = (unsigned char)(filesize>>24);

    bmpinfoheader[ 4] = (unsigned char)(       w    );
    bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
    bmpinfoheader[ 6] = (unsigned char)(       w>>16);
    bmpinfoheader[ 7] = (unsigned char)(       w>>24);
    bmpinfoheader[ 8] = (unsigned char)(       h    );
    bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
    bmpinfoheader[10] = (unsigned char)(       h>>16);
    bmpinfoheader[11] = (unsigned char)(       h>>24);

    FILE *f;
    f = fopen("temp.bmp","wb");
    fwrite(bmpfileheader,1,14,f);
    fwrite(bmpinfoheader,1,40,f);
    for (int i=0; i<h; i++) {
        long long offset = ((long long)w*(h-i-1)*3);
        fwrite(img+offset,3,w,f);
        fwrite(bmppad,1,(4-(w*3)%4)%4,f);
    }
    fclose(f);
}

unsigned char *createImage(State state) {
    int w = state.w;
    int h = state.h;

    if (w > MAX_WIDTH_HEIGHT) w = MAX_WIDTH_HEIGHT;
    if (h > MAX_WIDTH_HEIGHT) h = MAX_WIDTH_HEIGHT;

    if(rank!=0)
    {
      unsigned char r, g, b;
      int xstart = CHUNK_SIZE* (rank-1);
      unsigned char *img = NULL;
      if(img) free(img);
      long long size = (long long)w*(long long)CHUNK_SIZE*3;
      img = (unsigned char *) malloc(size);


      for(int px= xstart; px<xstart +CHUNK_SIZE;px++){
        for(int py = 0;py<h;py++){
                double xval = (px-w/2)/state.zoom + state.centerX;
                double yval = (py-h/2)/state.zoom + state.centerY;
                r = g = b = 0;
                float iterations = iterationsToEscape(xval,yval, state.maxIterations);
                if(iterations != -1){
                        float h = HUE_PER_ITERATION *iterations;
                        r= hue2rgb(h+120);
                        g= hue2rgb(h);
                        b= hue2rgb(h+240);
                }
                long long loc = (((long long)px-(long long)xstart)*(long long)w +(long long)py)*3;
                img[loc+2] = (unsigned char)(r);
                img[loc+1] = (unsigned char)(g);
                img[loc+0] = (unsigned char)(b);
        }
        }
        MPI_Send(&img[0],size,MPI_UNSIGNED_CHAR,0,0,MPI_COMM_WORLD);
        while(true){
                MPI_Recv(&xstart,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
                if(xstart==-1) break;
                for(int px= xstart; px<xstart +CHUNK_SIZE;px++){
        for(int py = 0;py<h;py++){
 double xval = (px-w/2)/state.zoom + state.centerX;
                double yval = (py-h/2)/state.zoom + state.centerY;
                r = g = b = 0;
                float iterations = iterationsToEscape(xval,yval, state.maxIterations);
                if(iterations != -1){
                        float h = HUE_PER_ITERATION *iterations;
                        r= hue2rgb(h+120);
                        g= hue2rgb(h);
                        b= hue2rgb(h+240);
                }
                long long loc = (((long long)px-(long long)xstart)*(long long)w +(long long)py)*3;
                img[loc+2] = (unsigned char)(r);
                img[loc+1] = (unsigned char)(g);
                img[loc+0] = (unsigned char)(b);
        }
        }
        MPI_Send(&img[0],size,MPI_UNSIGNED_CHAR,0,0,MPI_COMM_WORLD);
        }


    }
    if(rank==0)
    {
        unsigned char *imgReceive = NULL;
        if(imgReceive) free(imgReceive);
        long long size = (long long)w*(long long)CHUNK_SIZE*3;
        imgReceive = (unsigned char *) malloc(size);

        unsigned char *img = NULL;
        if(img) free(img);
        long long imgsize = (long long)w*(long long)h*3;
        img = (unsigned char *) malloc(imgsize);

        int startValues[PROCNUM];
        for(int i=0; i<PROCNUM; i++)
        {
                if(i!=0){
                        startValues[i] = CHUNK_SIZE*(i-1);
                }
        }

        int startingPoint = PROCNUM-1;
        for(int px=startingPoint;px<w;px+=CHUNK_SIZE){
                        int py = 0;
                        MPI_Recv(&imgReceive[0],size,MPI_UNSIGNED_CHAR,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
                        int receivedFrom = status.MPI_SOURCE;
                        int startingIndex = startValues[receivedFrom];
                        for(int rpx=startingIndex;rpx<startingIndex+CHUNK_SIZE;rpx++){
                                for(int rpy=0;rpy<h;rpy++){
                                        long long loc = ((long long)rpx+(long long)rpy*(long long)w)*3;
                                        long long rloc = (((long long)rpx-(long long)startingIndex)*(long long)w +(long long)py)*3;
                                        img[loc+2] = imgReceive[rloc+2];
                                        img[loc+1] = imgReceive[rloc+1];
                                        img[loc+0] = imgReceive[rloc+0];
                                }
                        }
                        MPI_Send(&px,1,MPI_INT,receivedFrom,0,MPI_COMM_WORLD);
                        startValues[receivedFrom] = px;
        }
        for(int i=0; i<PROCNUM-1;i++){
                int py = 0;
                        MPI_Recv(&imgReceive[0],size,MPI_UNSIGNED_CHAR,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
                        int receivedFrom = status.MPI_SOURCE;
                        int startingIndex = startValues[receivedFrom];
                        for(int rpx=startingIndex;rpx<startingIndex+CHUNK_SIZE;rpx++){
                                for(int rpy=0;rpy<h;rpy++){
                                        long long loc = ((long long)rpx+(long long)rpy*(long long)w)*3;
                                        long long rloc = (((long long)rpx-(long long)startingIndex)*(long long)w +(long long)py)*3;
                                        img[loc+2] = imgReceive[rloc+2];
   img[loc+1] = imgReceive[rloc+1];
                                        img[loc+0] = imgReceive[rloc+0];
                                }
                        }
                        int finish = -1;
                        MPI_Send(&finish,1,MPI_INT,receivedFrom,0,MPI_COMM_WORLD);
        }
        return img;
    }
    unsigned char *img = NULL;
    return img;
}
void draw(State state) {
    double start = When();
    unsigned char *img = createImage(state);
    double end = When();
    printf("It took %f seconds\n",end-start);
    if(rank==0)writeImage(img, state.w, state.h);
}

int main(int argc, char **argv) {

    //MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //printf("hello world from proc %d\n",rank);
    State state;
    draw(state);
    MPI_Finalize();
}

double When()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}



