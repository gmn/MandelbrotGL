/**
 *
 * Simple Mandelbrot Set viewer - (c) 2018, Gregory Naughton
 *
 */

// c headers
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#if defined(_WIN32) || defined(_WIN64)
 #include <windows.h>
#else
 #include <dirent.h>
#endif

// cpp headers
#include <string>
#include <iostream>

// OpenGL / glew Headers
#define GL3_PROTOTYPES 1
#include <GL/glew.h>
//#include <GL/gl.h>

// SDL2 Headers
//#include <SDL2/SDL.h>
#include <SDL.h>

// stb headers
#define STB_TRUETYPE_IMPLEMENTATION  // force following include to generate implementation
#include "stb_truetype.h"
#define FONT_FILENAME "../DejaVuSans-Bold.ttf"

// imgui headers
#include "imgui/imgui.h"
#include "imgui/imgui_impl_sdl.h"
#include "imgui/imgui_impl_opengl3.h"

#define GL_VERSION_MAJOR 3
#define GL_VERSION_MINOR 3

#define DEFAULT_WINDOW_WIDTH 1280
#define DEFAULT_WINDOW_HEIGHT 720


const unsigned int filename_zeros = 5;

static void usage( const char * exename ) {
    fprintf( stdout, "usage: %s [-xy WxH][-f]\n", exename );
    fprintf( stdout, "  -xy <width>x<height>    set the window resolution (higher looks better)\n" );
    fprintf( stdout, "  -f                      activate full screen\n" );
    fprintf( stdout, "example:\n  %s -f -xy 1920x1080\n", exename );
    fflush( stdout );
}

static bool _getint( const char * p, long int *li ) {
    if ( !p || !*p || !li ) {
        return false;
    }
    long int ival = strtol( p, 0, 10 );
    if ( errno == EINVAL || errno == ERANGE || ival == LONG_MAX || ival == LONG_MIN ) {
        return false;
    }
    *li = ival;
    return true;
}


#if defined(_WIN32) || defined(_WIN64)
int GetNextLowestFilenameNumber( const char * prefix ) {
	HANDLE hFind;
	WIN32_FIND_DATA FindFileData;
	int count = 0;
	size_t plen = strlen(prefix);
	const char * c;

	if ( (hFind = FindFirstFile( "./*.bmp", &FindFileData )) != INVALID_HANDLE_VALUE ) {
		do {
			printf( "Found: %s\n", FindFileData.cFileName );
			if ( (c = strstr(FindFileData.cFileName, prefix)) ) {
				c += plen;
				char buf[10];
				memcpy( buf, c, filename_zeros );
				buf[filename_zeros] = 0;
				long int x;
				if ( !_getint(buf, &x) )
					continue;
				if ( x > count )
					count = x;
			}
		} while ( FindNextFile( hFind, &FindFileData ) );
		FindClose( hFind );
	}
	return count + 1;
}

#else /* Not Windows */

// We write numbered filenames. We have to check for any with the current prefix so we dont over write them.
int GetNextLowestFilenameNumber( const char * prefix ) {
    DIR * dir_p = opendir( "." );
    if ( !dir_p ) {
        fprintf( stderr, "can't read directory: \".\"\n" );
        return -1;
    }

    struct dirent * ent = 0;
    int count = 0; // 0000
    size_t plen = strlen( prefix );
    const char * c;

    // count total files in directory, so that we can allocate
    while ( (ent = readdir( dir_p )) )
    {
        if ( *ent->d_name != '.' && ent->d_type == DT_REG ) {
            if ( (c = strstr( ent->d_name, prefix )) ) {
                c += plen;
                char buf[10];
                memcpy( buf, c, filename_zeros );
                buf[ filename_zeros ] = 0;
                long int x;
                if ( !_getint( buf, &x ) )
                    continue;
                if ( x > count )
                    count = x;
            }
        }
    }
    closedir( dir_p );

    return count + 1;
}
#endif

// IMGui helper
static void ShowHelpMarker(const char* desc)
{
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered())
    {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}


// writes one kind of Microsoft Bitmap file
struct BMPWriter {
#pragma pack( push, 1 )
    struct bmpFileHeader {
        char bfType[2];                 //  2   The characters "BM"
        unsigned int bfSize;            //  4   The size of the file in bytes
        unsigned short bfReserved1;     //  2   Unused - must be zero
        unsigned short bfReserved2;     //  2   Unused - must be zero
        unsigned int bfOffBits;         //  4   Offset to start of Pixel Data
#ifdef __GNUC__
    } __attribute__((packed));
#else
    };
#endif

    struct bmpImageHeader {
        unsigned int biSize;            //  4   Header Size - Must be at least 40
        unsigned int biWidth;           //  4   Image width in pixels
        unsigned int biHeight;          //  4   Image height in pixels
        unsigned short biPlanes;        //  2   Must be 1
        unsigned short biBitCount;      //  2   Bits per pixel - 1, 4, 8, 16, 24, or 32
        unsigned int biCompression;     //  4   Compression type (0 = uncompressed)
        unsigned int biSizeImage;       //  4   Image Size - may be zero for uncompressed images
        unsigned int biXPelsPerMeter;   //  4   Preferred resolution in pixels per meter
        unsigned int biYPelsPerMeter;   //  4   Preferred resolution in pixels per meter
        unsigned int biClrUsed;         //  4   Number Color Map entries that are actually used
        unsigned int biClrImportant;    //  4   Number of significant colors
#ifdef __GNUC__
    } __attribute__((packed));
#else
    };
#endif

    struct bmpNT4HeaderExtension {
        unsigned int RedMask;       /* Mask identifying bits of red component */
        unsigned int GreenMask;     /* Mask identifying bits of green component */
        unsigned int BlueMask;      /* Mask identifying bits of blue component */
        unsigned int AlphaMask;     /* Mask identifying bits of alpha component */
        unsigned int CSType;        /* Color space type */
        unsigned int RedX;          /* X coordinate of red endpoint */
        unsigned int RedY;          /* Y coordinate of red endpoint */
        unsigned int RedZ;          /* Z coordinate of red endpoint */
        unsigned int GreenX;        /* X coordinate of green endpoint */
        unsigned int GreenY;        /* Y coordinate of green endpoint */
        unsigned int GreenZ;        /* Z coordinate of green endpoint */
        unsigned int BlueX;         /* X coordinate of blue endpoint */
        unsigned int BlueY;         /* Y coordinate of blue endpoint */
        unsigned int BlueZ;         /* Z coordinate of blue endpoint */
        unsigned int GammaRed;      /* Gamma red coordinate scale value */
        unsigned int GammaGreen;    /* Gamma green coordinate scale value */
        unsigned int GammaBlue;     /* Gamma blue coordinate scale value */
#ifdef __GNUC__
    } __attribute__((packed));
#else
    };
#endif

#pragma pack(pop)

    static bool WriteRGB( const char * filename, unsigned char * pixels, int width, int height ) {
        bmpFileHeader header;
        bmpImageHeader image;
        size_t szFh = sizeof(bmpFileHeader);
        size_t szIh = sizeof(bmpImageHeader);
        size_t bufferSz = width * height * 3;

        memset( &header, 0, sizeof( bmpFileHeader ) );
        header.bfType[0] = 'B';
        header.bfType[1] = 'M';
        header.bfSize = szFh + szIh + bufferSz;
        header.bfOffBits = ( szFh + szIh );

        memset( &image, 0, sizeof( bmpImageHeader ) );
        image.biSize = szIh;
        image.biWidth = width;
        image.biHeight = height;
        image.biPlanes = 1;
        image.biBitCount = 24;
        image.biSizeImage = width * height * 3;
        image.biXPelsPerMeter = 0xb12;
        image.biYPelsPerMeter = 0xb12;
        image.biClrUsed = 1 << image.biBitCount;
        image.biClrImportant = 1 << image.biBitCount;

        FILE * fh = fopen( filename, "wb" );
        if ( !fh ) {
            return false;
        }

        fwrite( (void*) &header, szFh, 1, fh );
        fwrite( (void*) &image, szIh, 1, fh );

        int line = height - 1;
        line = 0;
        const int row_pixels = width * 3;
        while ( line < height ) {
            unsigned char * p = &pixels[ line * width * 3 ];
            for ( int i = 0; i < row_pixels; i+=3 ) {
                // flip R<-->B in RGB
                unsigned char t[3] = { p[ i + 2 ], p[ i + 1 ], p[ i ] };
                fwrite( (void*) t, 1, 3, fh );
            }
            ++line;
        }

        fclose(fh);
        return true;
    }
}; // struct BMPWriter


// naked base class, kind of like an interface
class Mandelbrot {
protected:
    int width;
    int height;

public:
    Mandelbrot( int w, int h ) : width(w), height(h) {
    }

    void Report( int start, int fin ) {
        Uint32 delta_ms = fin - start;
        unsigned int sec = delta_ms / 1000u;
        unsigned int tenths = ( delta_ms - sec ) / 100u;
        fprintf(stderr, "Mandelbrot Set FINISHED in %d.%d sec\n", sec, tenths );
        fflush(stderr);
    }

    virtual void Compute( unsigned char * pixel_buffer )
    {
        Uint32 start_ms = SDL_GetTicks();

        Uint32 finished_ms = SDL_GetTicks();
        Report( start_ms, finished_ms );
    }
}; // class Mandelbrot

// impl
class MandelbrotV2 : public Mandelbrot {
public:
    enum filter_t {
        FILTER_RGB,     // RGB == unchanged
        FILTER_RBG,     // the rest are channel swapped
        FILTER_GRB,
        FILTER_GBR,
        FILTER_BRG,
        FILTER_BGR
    };
    static const int TOTAL_COLOR_FILTERS = 6;

    struct IntDoubleDouble {
        int n;
        double Tr;
        double Ti;

        IntDoubleDouble() { }
        IntDoubleDouble( int a, double b, double c ) : n(a), Tr(b), Ti(c) { }
        void set( int a, double b, double c ) { n = a; Tr = b; Ti = c; }
    };

    struct color_t {
        unsigned char r, g, b;

        unsigned char operator[]( int i ) {
            return ( &r )[ i ];
        }
        void set( unsigned char r_, unsigned char g_, unsigned char b_ ) {
            r = r_; g = g_; b = b_;
        }
    };

private:

    double zoomStart;
    double lookAt[2];
    double zoom[2];

    // originally 360.0, this number determines how quickly the colors shift from one end of the spectrum to the other.
    double hue;
    double saturation;
    double vibrance; // blows out the gradient

    filter_t colorFilter;
    unsigned int ms_to_compute;

    void (MandelbrotV2::*colorPicker)( int, int, double, double, color_t * );

public:

    MandelbrotV2( int w, int h ) :
            Mandelbrot( w, h ),
            zoomStart( 3.14159265358979 ),
            lookAt{ -0.5, 0.0 },
            hue( 150.0 ),
            saturation( 1.0 ),
            vibrance( 14.0 ),
            colorFilter( FILTER_BGR ),
            ms_to_compute( 0 ),
            colorPicker( &MandelbrotV2::colorPicker_BGR )
    {
        zoom[0] = zoomStart;
        zoom[1] = zoomStart;
    }

    void SetDimensions( int w, int h ) {
        width = w;
        height = h;
    }

    void ResetValues() {
        lookAt[ 0 ] = -0.5;
        lookAt[ 1 ] = 0.0;
        zoom[ 0 ] = zoomStart;
        zoom[ 1 ] = zoomStart;
        hue = 150.0;
        saturation = 1.0;
        vibrance = 14.0;
        colorFilter = FILTER_BGR;
    }

    void ResetLocation() {
        lookAt[ 0 ] = -0.5;
        lookAt[ 1 ] = 0.0;
        zoom[ 0 ] = zoomStart;
        zoom[ 1 ] = zoomStart;
    }

    unsigned int GetComputeTime() { return ms_to_compute; }

    double * GetZoom() {
        return zoom;
    }

    double GetZoomStart() {
        return zoomStart;
    }

    void SetZoom( double new_zoom ) {
        zoom[0] = new_zoom;
        zoom[1] = new_zoom;
    }

    void ZoomInc( double zInc ) {
        zoom[0] += zInc;
        zoom[1] += zInc;
    }

    void ZoomMult( double zMult ) {
        zoom[0] *= zMult;
        zoom[1] *= zMult;
    }

    double * GetLookAt() {
        return lookAt;
    }

    void IncLookAt( double x, double y ) {
        lookAt[0] += x;
        lookAt[1] += y;
    }

    void SetLookAt( double x, double y ) {
        lookAt[0] = x;
        lookAt[1] = y;
    }

    void GetRange( double * xRange, double * yRange ) {
        xRange[0] = lookAt[0] - zoom[0] / 2.0;
        xRange[1] = lookAt[0] + zoom[0] / 2.0;
        yRange[0] = lookAt[1] - zoom[1] / 2.0;
        yRange[1] = lookAt[1] + zoom[1] / 2.0;
        AdjustAspectRatio( xRange, yRange );
    }

    void AdjustAspectRatio( double * xRange, double * yRange ) {
        double sratio = ((double)width)/(double)height;
        double ratio = fabs(xRange[1]-xRange[0]) / fabs(yRange[1]-yRange[0]);
        if ( sratio > ratio ) {
            xRange[0] *= ( sratio / ratio );
            xRange[1] *= ( sratio / ratio );
            zoom[0]   *= ( sratio / ratio );
        } else {
            yRange[0] *= ( ratio / sratio );
            yRange[1] *= ( ratio / sratio );
            zoom[1]   *= ( ratio / sratio );
        }
    }

    IntDoubleDouble * IterateEquation( double Cr, double Ci, double escapeRadius, int iterations, IntDoubleDouble * idd_p )
    {
        double Zr = 0.0;
        double Zi = 0.0;
        double Tr = 0.0;
        double Ti = 0.0;
        int    n  = 0;

        for ( ; n < iterations && (Tr + Ti) <= escapeRadius; ++n ) {
            Zi = 2 * Zr * Zi + Ci;
            Zr = Tr - Ti + Cr;
            Tr = Zr * Zr;
            Ti = Zi * Zi;
        }

        /*
         * Four more iterations (unrolled) to decrease error term;
         */
        //1
        Zi = 2 * Zr * Zi + Ci;
        Zr = Tr - Ti + Cr;
        Tr = Zr * Zr;
        Ti = Zi * Zi;

        //2
        Zi = 2 * Zr * Zi + Ci;
        Zr = Tr - Ti + Cr;
        Tr = Zr * Zr;
        Ti = Zi * Zi;

        //3
        Zi = 2 * Zr * Zi + Ci;
        Zr = Tr - Ti + Cr;
        Tr = Zr * Zr;
        Ti = Zi * Zi;

        //4
        Zi = 2 * Zr * Zi + Ci;
        Zr = Tr - Ti + Cr;
        Tr = Zr * Zr;
        Ti = Zi * Zi;

        //idd_p->set( n, Zr, Zi );
        idd_p->set( n, Tr, Ti );
        return idd_p;
    }

    // returns new value
    double HueUp() {
        double amt = 10.0;
        hue += amt;
        return hue;
    }

    // returns new value
    double HueDown() {
        double amt = 10.0;
        if ( hue > amt ) {
            hue -= amt;
        }
        return hue;
    }

    // returns new value
    double SaturationUp() {
        double amt = 0.05;
        saturation += amt;
        return saturation;
    }

    // returns new value
    double SaturationDown() {
        double amt = 0.05;
        if ( saturation > amt ) {
            saturation -= amt;
        }
        return saturation;
    }

    // returns new value
    double VibranceUp() {
        double amt = 0.5;
        vibrance += amt;
        return vibrance;
    }

    // returns new value
    double VibranceDown() {
        double amt = 0.5;
        if ( vibrance > amt ) {
            vibrance -= amt;
        }
        return vibrance;
    }

    double * GetHSV() {
        static double m[ 3 ];
        m[ 0 ] = hue;
        m[ 1 ] = saturation;
        m[ 2 ] = vibrance;
        return m;
    }

    void GetHue( double ** hue_ ) { *hue_ = &hue; }
    void GetSaturation( double ** saturation_ ) { *saturation_ = &saturation; }
    void GetVibrance( double ** vibrance_ ) { *vibrance_ = &vibrance; }

    // Some constants used with smoothColor
    const double logBase = 1.0 / log(2.0);
    const double logHalfBase = log(0.5) * logBase;

#if 0
    double smoothColor( int steps, double Tr, double Ti )
    {
        //return 1 + steps - log(log(sqrt(Zr*Zr+Zi*Zi)))/log(2.0);
        //return 1 + steps - log(log(sqrt(Tr*Tr+Ti*Ti)))/log(2.0);
        return 2.312156 + ((double)steps) - log(log(sqrt(Tr*Tr+Ti*Ti))) * logBase;
        //return 5.0 + ((double)steps) - logHalfBase - log( log(Tr + Ti) ) * logBase;
        //return 5 + steps - logHalfBase - log( log(Tr + Ti) ) * logBase;
    }
#endif

#define SMOOTH_COLOR( steps_, Tr_, Ti_ ) ( 2.312156 + ((double)steps_) - log(log(sqrt(Tr_*Tr_+Ti_*Ti_))) * logBase )

//#define SMOOTH_COLOR( steps_, Tr_, Ti_ ) ( 5.602512 + ((double)steps_) - logHalfBase - log( log(Tr_ + Ti_) ) * logBase )

/*
 * Convert hue-saturation-value/luminosity to RGB.
 *
 * Input ranges:
 *   H =   [0, 360] (integer degrees)
 *   S = [0.0, 1.0] (float)
 *   V = [0.0, 1.0] (float)
 */
#define HSV_TO_RGB( h__, s__, v__, color__ ) do {                                   \
    double h_ = h__;                                                                \
    double s_ = s__;                                                                \
    double v_ = v__;                                                                \
                                                                                    \
    if ( v_ > 1.0 )                                                                 \
        v_ = 1.0;                                                                   \
                                                                                    \
    double hp_ = h_ / 60.0;                                                         \
    double c_ = v_ * s_;                                                            \
                                                                                    \
    double rem_ = hp_ - floor(hp_);                                                 \
                                                                                    \
    double x_ = c_ * (1.0 - fabs( (((int)hp_) % 2) + rem_ - 1.0 ) );                \
    double rgb_[3] = {0.0, 0.0, 0.0};                                               \
                                                                                    \
    if ( 0.0 <= hp_ && hp_ < 1.0 ) { rgb_[0] = c_; rgb_[1] = x_; rgb_[2] = 0.0; }   \
    if ( 1.0 <= hp_ && hp_ < 2.0 ) { rgb_[0] = x_; rgb_[1] = c_; rgb_[2] = 0.0; }   \
    if ( 2.0 <= hp_ && hp_ < 3.0 ) { rgb_[0] = 0.0; rgb_[1] = c_; rgb_[2] = x_; }   \
    if ( 3.0 <= hp_ && hp_ < 4.0 ) { rgb_[0] = 0.0; rgb_[1] = x_; rgb_[2] = c_; }   \
    if ( 4.0 <= hp_ && hp_ < 5.0 ) { rgb_[0] = x_; rgb_[1] = 0.0; rgb_[2] = c_; }   \
    if ( 5.0 <= hp_ && hp_ < 6.0 ) { rgb_[0] = c_; rgb_[1] = 0.0; rgb_[2] = x_; }   \
                                                                                    \
    double m_ = v_ - c_;                                                            \
    rgb_[0] += m_;                                                                  \
    rgb_[1] += m_;                                                                  \
    rgb_[2] += m_;                                                                  \
                                                                                    \
    rgb_[0] *= 255.0;                                                               \
    rgb_[1] *= 255.0;                                                               \
    rgb_[2] *= 255.0;                                                               \
                                                                                    \
    color__->r = rgb_[0] > 255.0 ? 255 : rgb_[0];                                   \
    color__->g = rgb_[1] > 255.0 ? 255 : rgb_[1];                                   \
    color__->b = rgb_[2] > 255.0 ? 255 : rgb_[2];                                   \
} while(0)

#define SWAP_RED_AND_BLUE( color_ ) do {    \
            double t_ = color_->r;          \
            color_->r = color_->b;          \
            color_->b = t_;                 \
} while( 0 )

#define SWAP_GREEN_AND_BLUE( color_ ) do {  \
            double t_ = color_->g;          \
            color_->g = color_->b;          \
            color_->b = t_;                 \
} while( 0 )

#define SWAP_RED_AND_GREEN( color_ ) do {   \
            double t_ = color_->r;          \
            color_->r = color_->g;          \
            color_->g = t_;                 \
} while( 0 )


    void colorPicker_RGB( int max, int steps, double Tr, double Ti, color_t * color ) {
        if ( max == steps ) { // converged?
            color->set(0,0,0);
            return;
        }

        double v = SMOOTH_COLOR( steps, Tr, Ti );

        HSV_TO_RGB( (hue * v / ((double)max)), saturation, (vibrance * v / ((double)max)), color );
    }

    void colorPicker_RBG( int max, int steps, double Tr, double Ti, color_t * color ) {
        if ( max == steps ) { // converged?
            color->set(0,0,0);
            return;
        }

        double v = SMOOTH_COLOR( steps, Tr, Ti );

        HSV_TO_RGB( (hue * v / ((double)max)), saturation, (vibrance * v / ((double)max)), color );

        // R G B --> R B G
        SWAP_GREEN_AND_BLUE( color );
    }

    void colorPicker_GRB( int max, int steps, double Tr, double Ti, color_t * color ) {
        if ( max == steps ) { // converged?
            color->set(0,0,0);
            return;
        }

        double v = SMOOTH_COLOR( steps, Tr, Ti );

        HSV_TO_RGB( (hue * v / ((double)max)), saturation, (vibrance * v / ((double)max)), color );

        // R G B --> G R B
        SWAP_RED_AND_GREEN( color );
    }

    void colorPicker_GBR( int max, int steps, double Tr, double Ti, color_t * color ) {
        if ( max == steps ) { // converged?
            color->set(0,0,0);
            return;
        }

        double v = SMOOTH_COLOR( steps, Tr, Ti );

        HSV_TO_RGB( (hue * v / ((double)max)), saturation, (vibrance * v / ((double)max)), color );

        // R G B --> G R B
        SWAP_RED_AND_GREEN( color );

        // G R B --> G B R
        SWAP_GREEN_AND_BLUE( color );
    }

    void colorPicker_BRG( int max, int steps, double Tr, double Ti, color_t * color ) {
        if ( max == steps ) { // converged?
            color->set(0,0,0);
            return;
        }

        double v = SMOOTH_COLOR( steps, Tr, Ti );

        HSV_TO_RGB( (hue * v / ((double)max)), saturation, (vibrance * v / ((double)max)), color );

        // R G B --> G R B
        SWAP_RED_AND_GREEN( color );

        // G R B --> B R G
        SWAP_RED_AND_BLUE( color );
    }

    void colorPicker_BGR( int max, int steps, double Tr, double Ti, color_t * color ) {
        if ( max == steps ) { // converged?
            color->set(0,0,0);
            return;
        }

        double v = SMOOTH_COLOR( steps, Tr, Ti );

        HSV_TO_RGB( (hue * v / ((double)max)), saturation, (vibrance * v / ((double)max)), color );

        // R G B --> B G R
        SWAP_RED_AND_BLUE( color );
    }

    void colorPicker_greyScale( int max, int steps, double Tr, double Ti, color_t * c ) {
        if ( max == steps ) { // converged?
            c->set(0,0,0);
            return;
        }

        double dv = SMOOTH_COLOR( steps, Tr, Ti );
        int v = floor( (512.0 * dv) / ((double)steps) );
        if ( v > 255 )
            v = 255;
        else if ( v < 0 )
            v = 0;
        c->set( v, v, v );
    }

    void SetColorFilter( filter_t newFilter ) {
        if ( newFilter == colorFilter )
            return;
        colorFilter = newFilter;
        switch( newFilter ) {
        case FILTER_RGB:
            colorPicker = &MandelbrotV2::colorPicker_RGB;
            break;
        case FILTER_RBG:
            colorPicker = &MandelbrotV2::colorPicker_RBG;
            break;
        case FILTER_GRB:
            colorPicker = &MandelbrotV2::colorPicker_GRB;
            break;
        case FILTER_GBR:
            colorPicker = &MandelbrotV2::colorPicker_GBR;
            break;
        case FILTER_BRG:
            colorPicker = &MandelbrotV2::colorPicker_BRG;
            break;
        case FILTER_BGR:
            colorPicker = &MandelbrotV2::colorPicker_BGR;
            break;
        }
    }

    filter_t GetColorFilter() {
        return colorFilter;
    }

    void PrevFilter() {
        int x = (int) colorFilter;
        if ( --x < 0 )
            x = TOTAL_COLOR_FILTERS - 1;
        SetColorFilter( (filter_t) x );
    }

    void NextFilter() {
        int x = (int) colorFilter;
        if ( ++x >= TOTAL_COLOR_FILTERS )
            x = 0;
        SetColorFilter( (filter_t) x );
    }

    virtual void Compute( unsigned char * pixel_buffer ) {
        Uint32 start_ms = SDL_GetTicks();

        /*
        double xRange[2] = { lookAt[0]-zoom[0]/2.0, lookAt[0]+zoom[0]/2.0 };
        double yRange[2] = { lookAt[1]-zoom[1]/2.0, lookAt[1]+zoom[1]/2.0 };
        AdjustAspectRatio( xRange, yRange );
        */
        double xRange[2];
        double yRange[2];
        GetRange( xRange, yRange );

        double escapeRadius = 10 * 10;

//        double dx = ( xRange[1] - xRange[0] ) / ((double) width);
//        double dy = ( yRange[1] - yRange[0] ) / ((double)height);
        double dx = ( xRange[1] - xRange[0] ) / ( 0.5 + (((double) width)-1.0) );
        double dy = ( yRange[1] - yRange[0] ) / ( 0.5 + (((double)height)-1.0) );

        double Cr = xRange[0];
        double Ci = yRange[0];

        int max_iterations = floor( 223.0 / sqrt( 0.001 + 2.0 * fmin( fabs( xRange[0] - xRange[1] ), fabs( yRange[0] - yRange[1] ) ) ) );

        int offset = 0;

        // for each line
        int sy = 0;
        while ( sy < height )
        {
            Cr = xRange[0];

            for ( int x = 0 ; x < width; ++x, Cr += dx ) {
                struct IntDoubleDouble idd;
                IterateEquation( Cr, Ci, escapeRadius, max_iterations, &idd );
                color_t c;
                (*this.*colorPicker)( max_iterations, idd.n, idd.Tr, idd.Ti, &c );
                pixel_buffer[ offset++ ] = c[0];
                pixel_buffer[ offset++ ] = c[1];
                pixel_buffer[ offset++ ] = c[2];
            }

            Ci += dy;
            sy++;
        }

        Uint32 finished_ms = SDL_GetTicks();

        ms_to_compute = finished_ms - start_ms;
        Report( start_ms, finished_ms );
    }

}; // class MandelbrotV2


class TimedMessage {
private:
    char buf[ 1024 ];
    unsigned int t0;

public:
    TimedMessage() : t0( 0 ) {
    }

    void Set( const char * msg, unsigned int ttl, int x_, int y_ ) {
        strncpy( buf, msg, sizeof(buf) );
        t0 = SDL_GetTicks() + ttl;
        x = x_;
        y = y_;
    }

    bool IsAlive() {
        unsigned int now = SDL_GetTicks();
        return ( now < t0 );
    }

    char * c_str() { return buf; }

public:
    int x; int y;
}; // class TimedMessage


// handles application level behavior, ie, window size, key-input, sdl & opengl setup.
// passes user input to mandelbrot calculator
class Application {
private:

    int windowWidth;
    int windowHeight;
    int bytesPerPixel;

    std::string windowTitle;

    // Our SDL_Window ( just like with SDL2 wihout OpenGL)
    SDL_Window * mainWindow;

    // Our opengl context handle
    SDL_GLContext mainGLContext;

    unsigned char * screenTextureBuffer;

    // for the persistent mandelbrot texture
    GLuint textureHandle;
    GLuint arrayHandle;

    bool drawInside;
    bool drawOutside;
    bool dataChanged;
    int drawCrosshair;

    bool fullscreen_flag;

    MandelbrotV2 mandelbrot;

    // STB font buffers
    stbtt_bakedchar * cdata = nullptr;
    GLuint * ftex = nullptr;
    int font_scale;
    float font_px;
    char * message_p;

    int * mouse_line;
    int mouse_active;
    int mouse_activity_start;

    TimedMessage fileMessage;

    int currentMouseLocation[2];
    struct desktopInfo_t {
        int width, height, refreshRate;
        int dummy;
    } desktopInfo;

    bool gui_show_control_panel;
    bool gui_show_help_window;
    bool gui_show_demo_window;
    bool gui_quit_requested;
    bool gui_redraw_requested;
    bool gui_reset_all_requested;
    bool gui_reset_location_requested;
    bool gui_screenshot_requested;

    bool useDesktopResolution;
    bool gl_filter_on;
    GLuint m_shaderColor;
    GLuint m_shaderTexture;

    int m_AttribTexColors;
    int m_AttribPosition ;
    int m_AttribTexCoord ;
    int m_AttribInVertex ;
    int m_AttribInColor  ;

public:

    explicit Application( int width_, int height_, bool fullscreen_ =false, bool useDesktopResolution_ =true ) :
            bytesPerPixel( 3 ),
            mainWindow( nullptr ),
            screenTextureBuffer( nullptr ),
            drawInside( true ),
            drawOutside( true ),
            dataChanged( true ),
            drawCrosshair( 3 ),
            fullscreen_flag( 0 ),
            mandelbrot( width_, height_ ),
            font_scale( 256 ),
            font_px( 14.5 ),
            message_p( nullptr ),
            mouse_line( nullptr ),
            mouse_active( 0 ),
            mouse_activity_start( 0 ),
            gui_show_control_panel( true ),
            gui_show_help_window( false ),
            gui_show_demo_window( false ),
            gui_quit_requested( false ),
            gui_redraw_requested( false ),
            gui_reset_all_requested( false ),
            gui_reset_location_requested( false ),
            gui_screenshot_requested( false ),
            useDesktopResolution( useDesktopResolution_ ),
            gl_filter_on( true )
    {
        windowTitle = "OpenGL Mandelbrot Viewer";

        windowWidth = width_;
        windowHeight = height_;
        if ( fullscreen_ ) {
            fullscreen_flag = SDL_WINDOW_FULLSCREEN;
        }
        currentMouseLocation[ 0 ] = width_ / 2;
        currentMouseLocation[ 1 ] = height_ / 2;
        memset( &desktopInfo, 0, sizeof( desktopInfo_t ) );
    }

    ~Application() {
        my_stbtt_cleanup();
        Shutdown();
    }

    // needs macro
    static void CheckSDLError(int line = -1)
    {
        std::string error = SDL_GetError();

        if (error != "")
        {
            std::cout << "SLD Error : " << error << std::endl;

            if (line != -1)
                std::cout << "\nLine : " << line << std::endl;

            SDL_ClearError();
        }
    }

    bool SetOpenGLAttributes()
    {
        // Set our OpenGL version.
        // 3.2 is part of the modern versions of OpenGL, but most video cards should be able to run it
        SDL_GL_SetAttribute( SDL_GL_CONTEXT_MAJOR_VERSION, GL_VERSION_MAJOR );
        SDL_GL_SetAttribute( SDL_GL_CONTEXT_MINOR_VERSION, GL_VERSION_MINOR );

        // Turn on double buffering with a 24bit Z buffer.
        // You may need to change this to 16 or 32 for your system
        SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
        SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 24 );
        SDL_GL_SetAttribute( SDL_GL_STENCIL_SIZE, 8 );

        return true;
    }

    static void GL_ResetOrthographicProjection( int width, int height )
    {
        if ( height == 0 )
            height = 1;

        glViewport( 0, 0, ( GLint )width, ( GLint )height );

        //
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity( );
        glOrtho( 0, width, 0, height, -1, 1 );

        //
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity ();
        glRasterPos3i( 0, 0, 0 );
        glTranslatef ( 0.375f, 0.375f, 0.f );
    }

    static int GL_InitState( void )
    {
        // specify color to clear to
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

        // specify the clear value of the depth buffer
        glClearDepth( 1.0f );
        glDepthFunc( GL_ALWAYS );
        glDisable( GL_DEPTH_TEST );
        //glEnable( GL_DEPTH_TEST );

        // blend
        //glDisable( GL_BLEND ); // blend incoming colors with colors in the buffers
        glEnable( GL_BLEND ); // blend incoming colors with colors in the buffers
        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        glBlendEquation( GL_FUNC_ADD );

        // cull backfaces
        glEnable(GL_CULL_FACE);

        // flat colors
        glShadeModel( GL_FLAT );

        glDisable( GL_ALPHA_TEST );
        glAlphaFunc( GL_LESS, 1.0 );

        return 1;
    }

    static void Print_SDL_GL_Attributes()
    {
        int value = 0;
        SDL_GL_GetAttribute( SDL_GL_CONTEXT_MAJOR_VERSION, &value );
        std::cout << "SDL_GL_CONTEXT_MAJOR_VERSION: " << value << std::endl;

        SDL_GL_GetAttribute( SDL_GL_CONTEXT_MINOR_VERSION, &value );
        std::cout << "SDL_GL_CONTEXT_MINOR_VERSION: " << value << std::endl;
    }

    void ImguiInit() {

        // Setup Dear ImGui binding
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO();
        (void)io;
        //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls

        ImGui_ImplSDL2_InitForOpenGL( mainWindow, mainGLContext );
        ImGui_ImplOpenGL3_Init();

        // Setup style
        //ImGui::StyleColorsLight();
        //ImGui::StyleColorsDark();
        ImGui::StyleColorsClassic();
    }

    bool Init()
    {
        // Initialize SDL's Video subsystem
        if ( SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) < 0 ) {
            std::cout << "Failed to init SDL\n";
            return false;
        }

        // do before window init
        SetOpenGLAttributes();

        // desktop display mode before running this app
        SDL_DisplayMode current;
        SDL_GetCurrentDisplayMode( 0, &current );
        desktopInfo.width = current.w;
        desktopInfo.height = current.h;
        desktopInfo.refreshRate = current.refresh_rate;

        // prefer current desktop resolution for fullscreen if user didn't specify one
        if ( fullscreen_flag && useDesktopResolution ) {
            windowWidth = desktopInfo.width;
            windowHeight = desktopInfo.height;
            currentMouseLocation[ 0 ] = desktopInfo.width / 2;
            currentMouseLocation[ 1 ] = desktopInfo.height / 2;
            mandelbrot.SetDimensions( desktopInfo.width, desktopInfo.height );
        }

        // create application
        if ( fullscreen_flag )
            printf( "Going fullscreen @ %d x %d.\n", windowWidth, windowHeight );
        else
            printf( "Creating %d x %d window.\n", windowWidth, windowHeight );

        // Create our window centered
        mainWindow = SDL_CreateWindow( windowTitle.c_str(), SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, windowWidth, windowHeight, SDL_WINDOW_OPENGL /*|SDL_WINDOW_RESIZABLE*/ );

        // Do we have a main window?
        if ( !mainWindow ) {
            std::cout << "Unable to create window\n";
            CheckSDLError(__LINE__); // FIXME I hate this naive function
            return false;
        }

        // Create our opengl context and attach it to our window
        mainGLContext = SDL_GL_CreateContext( mainWindow );

        // enable vsync
        SDL_GL_SetSwapInterval(1);

        // Init GLEW
        // Apparently, this is needed for Apple. Thanks to Ross Vander for letting me know
        #ifndef __APPLE__
        glewExperimental = GL_TRUE;
        glewInit();
        #endif
        //gl3wInit();

        ImguiInit();

        // allocate userspace texture buffer
        screenTextureBuffer = (unsigned char *) calloc( (unsigned int)( windowWidth * windowHeight * bytesPerPixel ) + 1u, 1u );

        InitScreenTexture( &textureHandle, &arrayHandle, windowWidth, windowHeight, GL_NEAREST );

        if ( !CreateShaderProgram_Identity( &m_shaderColor ) ) {
            return false;
        }
        if ( !CreateShaderProgram_Texture( &m_shaderTexture ) ) {
            return false;
        }

        SDL_ShowCursor(1);

        // FIXME: must redo for >= gl3.1
        GL_ResetOrthographicProjection( windowWidth, windowHeight );

        GL_InitState();

        my_stbtt_initfont();

        // fix scaling early - or else first click goes to wrong location
        ScalingFix();

        return true;
    }

    void ScalingFix() {
        CenterAtPoint( windowWidth/2, windowHeight/2 );
        double xRange[2];
        double yRange[2];
        mandelbrot.GetRange( xRange, yRange );
    }

    static void InitScreenTexture( GLuint * tex_obj, GLuint * array_handle, int width, int height, GLenum filterMode )
    {
        //
        // TEXTURE
        //
        glEnable( GL_TEXTURE_2D );

        // ask API state-machine to generate texture object handle
        glGenTextures( 1, tex_obj );

        // all subsequent commands affect this texture unit;
        // the active texture unit is inferred from this selector(enum)
        glActiveTexture( GL_TEXTURE0 );

        // set the generated id as a texture type
        glBindTexture( GL_TEXTURE_2D, *tex_obj );

        // verify it
        assert( glIsTexture( *tex_obj ) );

        //
        glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );

        // For immutable textures in opengl-3.3, use this with glTexSubImage2D()
        glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, (void*)0 );

        // GL_NEAREST displays the actual pixel from texture data rather than interpolation
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filterMode );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filterMode );

        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER );

        glBindTexture( GL_TEXTURE_2D, 0 );
        glDisable( GL_TEXTURE_2D );
    }

    bool CompileShaderPair( GLuint * progID, const char * vertSrc, const char * fragSrc )
    {
        // check compilation
        GLint status_return = 0;

        // create the shader units
        GLuint vertShader = glCreateShader( GL_VERTEX_SHADER );
        GLuint fragShader = glCreateShader( GL_FRAGMENT_SHADER );

        // associate each shader handle with its source code
        glShaderSource( vertShader, 1, (const GLchar **)&vertSrc, 0 );
        glShaderSource( fragShader, 1, (const GLchar **)&fragSrc, 0 );

        // compile
        glCompileShader( vertShader );
        glGetShaderiv( vertShader, GL_COMPILE_STATUS, &status_return );
        if ( GL_FALSE == status_return ) {
            fprintf( stderr, "failed to compile vertex shader\n" );
            return false;
        }

        glCompileShader( fragShader );
        glGetShaderiv( fragShader, GL_COMPILE_STATUS, &status_return );
        if ( GL_FALSE == status_return ) {
            fprintf( stderr, "failed to compile fragment shader\n" );
            return false;
        }

        // link together and create a gl shader program
        GLuint shaderProgram = glCreateProgram();

        glAttachShader( shaderProgram, vertShader );
        glAttachShader( shaderProgram, fragShader );

        glLinkProgram( shaderProgram );
        glGetProgramiv( shaderProgram, GL_LINK_STATUS, &status_return );
        if ( GL_FALSE == status_return ) {
            fprintf( stderr, "failed to link shaders\n" );
            return false;
        }

        *progID = shaderProgram;
        return true;
    }

    // passthru shaders
    bool CreateShaderProgram_Identity( GLuint * programID )
    {
        // identity vertex shader
        const char * vertex_shader =
            "#version 150\n"
            "in vec4 in_vertex;\n"
            "in vec4 in_color;\n"
            "out vec4 varying_color;\n"
            "void main(void) {\n"
            "    varying_color = in_color;\n"
            "    gl_Position = in_vertex;\n"
            "}\n";

        // identity fragment shader
        const char * fragment_shader =
            "#version 150\n"
            "// It was expressed that some drivers required this next line to function properly\n"
            "precision highp float;\n"
            "out vec4 out_fragColor;\n"
            "in  vec4 varying_color;\n"
            "void main(void)\n"
            "{\n"
            "    out_fragColor = varying_color;\n"
            "}\n";

        if ( !CompileShaderPair( programID, vertex_shader, fragment_shader ) ) {
            fprintf( stderr, "Failure in Identity Shader\n" );
            return false;
        }

        //glBindAttribLocation( *programID, 0, "in_vertex" );
        //glBindAttribLocation( *programID, 1, "in_color" );

        m_AttribInVertex  = glGetAttribLocation( *programID, "in_vertex" );
        m_AttribInColor  = glGetAttribLocation( *programID, "in_color" );

        return true;
    }

    bool CreateShaderProgram_Texture( GLuint * programID )
    {
        const char * tex_vert_shader =
            "#version 150\n"
            "in vec4 in_position;\n"
            "in vec2 in_texCoord;\n"
            "out vec2 inter_texCoord;\n"
            "void main(void)\n"
            "{\n"
            "   inter_texCoord = in_texCoord;\n"
            "   gl_Position = in_position;\n"
            "}\n";

        const char * tex_frag_shader =
            "#version 150\n"
            "precision highp float;\n"
            "out vec4 out_fragColor;\n"
            "uniform sampler2D texColors;\n"
            "in vec2 inter_texCoord;\n"
            "void main(void)\n"
            "{\n"
            "   out_fragColor = texture( texColors, inter_texCoord );\n"
            "}\n";

        //glBindAttribLocation( shaderProgram, 0, "in_position" );
        //glBindAttribLocation( shaderProgram, 1, "in_texCoord" );

        if ( !CompileShaderPair( programID, tex_vert_shader, tex_frag_shader ) ) {
            fprintf( stderr, "Failure in Texture Shader\n" );
            return false;
        }

        m_AttribTexColors = glGetUniformLocation( *programID, "texColors" );
        m_AttribPosition  = glGetAttribLocation( *programID, "in_position" );
        m_AttribTexCoord  = glGetAttribLocation( *programID, "in_texCoord" );

        return true;
    }


    void DrawScreenTexture( GLuint texNum, unsigned char * bytes )
    {
        glUseProgram( m_shaderTexture );

        glEnable( GL_TEXTURE_2D );

        //
        glActiveTexture( GL_TEXTURE0 );

        //load it into the graphics hardware:
        glBindTexture( GL_TEXTURE_2D, texNum );

        glUniform1i( m_AttribTexColors, 0 );
        if ( glBindSampler )
            glBindSampler( 0, 0 );

        //
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, gl_filter_on ? GL_LINEAR : GL_NEAREST );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, gl_filter_on ? GL_LINEAR : GL_NEAREST );

        // specify data for texture unit
        glTexSubImage2D( GL_TEXTURE_2D,
                         0,                                     // 1st mip-map level
                         0, 0,                                  // offset into texture
                         windowWidth, windowHeight,
                         GL_RGB, GL_UNSIGNED_BYTE, bytes );

        // specify texture vertices & and tex-coords through VBO
        static const GLfloat quad_data[] = {
            // verteces
            -1.0f, -1.0f, 0.0f, 1.0f,
            +1.0f, -1.0f, 0.0f, 1.0f,
            +1.0f, +1.0f, 0.0f, 1.0f,
            -1.0f, +1.0f, 0.0f, 1.0f,
            // tex coords
            0.0f, 0.0f,
            1.0f, 0.0f,
            1.0f, 1.0f,
            0.0f, 1.0f
        };

        // array buffer
        GLuint buf;
        glGenBuffers( 1, &buf );
        glBindBuffer( GL_ARRAY_BUFFER, buf );
        glBufferData( GL_ARRAY_BUFFER, sizeof(quad_data), quad_data,
                        GL_STATIC_DRAW /* modified once and used many times */
                    );

        // vertex array
        GLuint vao;
        glGenVertexArrays( 1, &vao );
        glBindVertexArray( vao );
        glVertexAttribPointer( 0,               // index of generic attribute(?)
                               4,               // number of components per vertex attribute
                               GL_FLOAT,        // component type
                               GL_FALSE,        // is normalized?
                               0,               // stride
                               (GLvoid*)0 );

        glEnableVertexAttribArray( m_AttribPosition );
        glVertexAttribPointer( 1, 2, GL_FLOAT, GL_FALSE, 0, (GLvoid*)( 16 * sizeof(GLfloat) ) );
        glEnableVertexAttribArray( m_AttribTexCoord );

        glDrawArrays( GL_TRIANGLE_FAN, 0, 4 );

        glBindTexture( GL_TEXTURE_2D, 0 );
        glDisable( GL_TEXTURE_2D );
    }

    void DrawCrosshairs()
    {
        static int lastCrosshair = drawCrosshair;
        int pixels = 10;

#define CROSSHAIR_MSG() do {\
    if ( drawCrosshair != lastCrosshair ) {\
        switch( drawCrosshair ) {\
        case 0: printf( "crosshair disabled\n" ); break;\
        case 1: printf( "small, white crosshair\n" ); break;\
        case 2: printf( "large, white crosshair\n" ); break;\
        case 3: printf( "small, grey crosshair\n" ); break;\
        case 4: printf( "large, grey crosshair\n" ); break;\
        }\
        lastCrosshair = drawCrosshair;\
    }\
}while(0)

        glUseProgram( m_shaderColor );

        GLfloat * color4p = nullptr;
        GLfloat white[] = {
            1.0f, 1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f, 1.0f
        };
        GLfloat grey[] = {
            0.4f, 0.4f, 0.4f, 1.0f,
            0.4f, 0.4f, 0.4f, 1.0f,
            0.4f, 0.4f, 0.4f, 1.0f,
            0.4f, 0.4f, 0.4f, 1.0f
        };

        if ( 0 == drawCrosshair ) {
            CROSSHAIR_MSG();
            return;
        } else if ( 1 == drawCrosshair ) {
            CROSSHAIR_MSG();
            pixels = 12;
            color4p = white;
        } else if ( 2 == drawCrosshair ) {
            CROSSHAIR_MSG();
            pixels = 32;
            color4p = white;
        } else if ( 3 == drawCrosshair ) {
            CROSSHAIR_MSG();
            pixels = 12;
            color4p = grey;
        } else if ( 4 == drawCrosshair ) {
            CROSSHAIR_MSG();
            pixels = 32;
            color4p = grey;
        }
#undef CROSSHAIR_MSG

        assert( color4p != nullptr );

        // use this when glOrtho is configured to screen resolution
        GLfloat w = pixels;
        GLfloat h = pixels;
        GLfloat c[2] = { (GLfloat)windowWidth/2.0f, (GLfloat)windowHeight/2.0f };

        // -1.0 <= x <= +1.0 ; -1.0 <= y <= +1.0
        c[0] = 0.0f;
        c[1] = 0.0f;
        w = pixels/(GLfloat)(windowWidth/2.0f);
        h = pixels/(GLfloat)(windowHeight/2.0f);

        const GLfloat line_cross[4][4] = {
            { c[0]-w, c[1], 0.0, 1.0 },
            { c[0]+w, c[1], 0.0, 1.0 },
            { c[0], c[1]-h, 0.0, 1.0 },
            { c[0], c[1]+h, 0.0, 1.0 }
        };

        glLineWidth( 1.0 );

        GLuint vao;
        GLuint vbo[2];

        glGenVertexArrays( 1, &vao );
        glBindVertexArray( vao );
        glGenBuffers( 2, vbo );

        // vertex array
        glBindBuffer( GL_ARRAY_BUFFER, vbo[0] );
        glBufferData( GL_ARRAY_BUFFER, sizeof(line_cross), line_cross, GL_STATIC_DRAW );
        glVertexAttribPointer( m_AttribInVertex, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0 );
        glEnableVertexAttribArray( m_AttribInVertex );

        // color array
        glBindBuffer( GL_ARRAY_BUFFER, vbo[1] );
        glBufferData( GL_ARRAY_BUFFER, sizeof(GLfloat)*16, color4p, GL_STATIC_DRAW );
        glVertexAttribPointer( m_AttribInColor, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0 );
        glEnableVertexAttribArray( m_AttribInColor );

        glDrawArrays( GL_LINES, 0, 4 );
    }

    void DrawMouseline() {
        if ( !mouse_line )
            return;

        // turn off mouse line if we're not active for N seconds
        int t1 = SDL_GetTicks();
        if ( 0 == mouse_activity_start || t1 - mouse_activity_start > 1500 ) {
            mouse_activity_start = t1;
            mouse_active = 0;
            mouse_line = nullptr;
            return;
        }

        glLineWidth( 1.0 );
        glBegin( GL_LINES );
            glVertex2i( windowWidth/2, windowHeight/2 );
            glVertex2i( mouse_line[0], mouse_line[1] );
        glEnd();
    }

    void ResetMouseline() {
        mouse_active = 0;
        mouse_line = nullptr;
        mouse_activity_start = 0;
    }

    void HandleMouseMove( int x, int y )
    {
        //printf( "x: %5d, y: %5d\n", x, y );
        static int xy[2];

        mouse_activity_start = SDL_GetTicks();

        xy[0] = x;
        xy[1] = windowHeight - y;

        if ( ++mouse_active > 6 ) {
            mouse_line = xy; // this instructs draw function to draw mouse line
        }

        currentMouseLocation[ 0 ] = x;
        currentMouseLocation[ 1 ] = y;
    }

    void HandleMouseWheel( int x, int y ) {
        printf( " --> mousewheel : (%d, %d)\n", x, y );
        //CenterAtPoint( currentMouseLocation[0], currentMouseLocation[1] );
        mandelbrot.ZoomMult( 1.0 + 0.142857 * y );
        dataChanged = true;
    }

    void SaveScreenshot() {
        char filename_prefix[ 256 ];
        sprintf( filename_prefix, "mandelbrot-%dx%d-", windowWidth, windowHeight );

        int img_counter = GetNextLowestFilenameNumber( filename_prefix );
        if ( img_counter < 0 ) {
            fprintf( stderr, "there was an error reading directory\n" );
            fflush(stderr);
            return;
        }

        char filename[256];
        snprintf( filename, sizeof(filename), "%s%05d.bmp", filename_prefix, img_counter );

        printf( "writing hardcopy to \"%s\"\n", filename );
        if ( !BMPWriter::WriteRGB( filename, screenTextureBuffer, windowWidth, windowHeight ) ) {
            fprintf( stderr, "there was an error writing to %s\n", filename );
            fflush(stderr);
            return;
        }

        char buf[1024];
        snprintf( buf, sizeof(buf), "wrote:  %s", filename );
        fileMessage.Set( buf, 3000, windowWidth/2 - 150, windowHeight/2 - 150 );
    }

    void CenterAtPoint( int x, int y ) {
        double * lookAt = mandelbrot.GetLookAt();
        double * zoom = mandelbrot.GetZoom();

        // mouse clicked this point in screen-space
        double point[2] = { (double)x, (double) y };
        double window[2] = { (double)windowWidth, (double)windowHeight };

        // conversion between window-space and mandelbrot-space
        double ratio[2] = { zoom[0] / window[0], zoom[1] / window[1] };

        // mouse click point in mandelbrot-space
        double p_scaled[2] = { ratio[0] * point[0], ratio[1] * point[1] };

        double xRange[2] = { lookAt[0]-zoom[0]/2.0, lookAt[0]+zoom[0]/2.0 };
        double yRange[2] = { lookAt[1]-zoom[1]/2.0, lookAt[1]+zoom[1]/2.0 };
        double center[2] = { (xRange[1]-xRange[0])/2.0 , (yRange[1]-yRange[0])/2.0 };

        // using vector arithmetic:
        //  center = p_scaled + q ; To shift to p we add -q ; -q = -1 * ( center - p_scaled )
        //mandelbrot.IncLookAt( -1.0 * ( center[0] - p_scaled[0] ), -1.0 * ( center[1] - p_scaled[1] ) );
        mandelbrot.IncLookAt( -1.0 * ( center[0] - p_scaled[0] ), ( center[1] - p_scaled[1] ) );
    }

    void HandleMouseClick( int x, int y, SDL_MouseButtonEvent * b_evt )
    {
      #define P(r) printf( " --> %s : (%d, %d)\n", #r, x, y );

        if ( b_evt->type == SDL_MOUSEBUTTONDOWN )
        {
            switch( b_evt->button )
            {
            case SDL_BUTTON_LEFT: {
                P(leftclick);

                ResetMouseline();

                CenterAtPoint( x, y );

                double * lookAt = mandelbrot.GetLookAt();
                printf( "lookAt: (%lf, %lf)\n", lookAt[0], lookAt[1] );
                dataChanged = true;
                break;
            }

            case SDL_BUTTON_MIDDLE:
                P(middleclick);
                break;
            case SDL_BUTTON_RIGHT:
                P(rightclick);
                break;
            }
        }
        else if ( b_evt->type == SDL_MOUSEBUTTONUP )
        {
        }

      #undef P
    }

    bool HandleKeyPress( SDL_Keysym * keysym )
    {
        switch ( keysym->sym )
        {
        case SDLK_ESCAPE:
        case SDLK_q:
            return false;
            break;
        case SDLK_f:
            fullscreen_flag = SDL_WINDOW_FULLSCREEN == fullscreen_flag ? 0 : SDL_WINDOW_FULLSCREEN;
            SDL_SetWindowFullscreen( mainWindow, fullscreen_flag );
            break;
        case SDLK_SPACE:
		case SDLK_RETURN:
            dataChanged = true;
            break;
        case SDLK_h:
            if ( ( KMOD_SHIFT & keysym->mod ) || ( KMOD_LSHIFT & keysym->mod ) ) {
                double d = mandelbrot.HueDown();
                printf( "hue down, %lf\n", d );
                dataChanged = true;
            } else {
                double d = mandelbrot.HueUp();
                printf( "hue up, %lf\n", d );
                dataChanged = true;
            }
            break;
        case SDLK_s:
            if ( ( KMOD_SHIFT & keysym->mod ) || ( KMOD_LSHIFT & keysym->mod ) ) {
                double d = mandelbrot.SaturationDown();
                printf( "saturation down, %lf\n", d );
                dataChanged = true;
            } else {
                double d = mandelbrot.SaturationUp();
                printf( "saturation up, %lf\n", d );
                dataChanged = true;
            }
            break;
        case SDLK_v:
            if ( ( KMOD_SHIFT & keysym->mod ) || ( KMOD_LSHIFT & keysym->mod ) ) {
                double d = mandelbrot.VibranceDown();
                printf( "vibrance down, %lf\n", d );
                dataChanged = true;
            } else {
                double d = mandelbrot.VibranceUp();
                printf( "vibrance up, %lf\n", d );
                dataChanged = true;
            }
            break;
        case SDLK_o:
            SaveScreenshot();
            break;

        case SDLK_m: {
            if ( ( KMOD_SHIFT & keysym->mod ) || ( KMOD_LSHIFT & keysym->mod ) ) {
                mandelbrot.PrevFilter();
            } else {
                mandelbrot.NextFilter();
            }
            dataChanged = true;
            break;
        }

        case SDLK_c:
            if ( 5 == ++drawCrosshair )
                drawCrosshair = 0;
            break;

        case SDLK_g:
            gui_show_control_panel = !gui_show_control_panel;
            break;

        case SDLK_QUESTION:
        case SDLK_SLASH:
            gui_show_help_window = !gui_show_help_window;
            break;

        case SDLK_w:
            gl_filter_on = !gl_filter_on;
            break;

        default:
            break;
        }

        return true;
    }

    // returns true if event was consumed by IMGUI
    bool ImguiProcessedEvent( SDL_Event * event_p ) {
        //bool mouseDragging = ImGui::GetCurrentContext()->MovingWindow != NULL;
        //bool mouseInsideGuiWindowBorder = false;

        ImGui_ImplSDL2_ProcessEvent( event_p );

        //return ( io.WantCaptureMouse || io.WantCaptureKeyboard );
        /**
         * Only block for keyboard input. This way 'q' key will still fall through and quit
         * when the mouse is hovering over Imgui. We only want to block keyboard when we
         * actually entering an input into a Scaler field.
         *
         * Also, we want the mouse to update even if its over the Gui window.
         */

        SDL_Keysym * key = &event_p->key.keysym;
        bool Quit = key->sym == SDLK_q || key->sym == SDLK_ESCAPE || event_p->type == SDL_QUIT;

        return !Quit && ( ImGui::GetIO().WantCaptureKeyboard || ImGui::GetIO().WantCaptureMouse );
    }

    void PollEvents( bool & loop, int & input_latency )
    {
        loop = true;

        SDL_Event event;
        while ( SDL_PollEvent( &event ) )
        {
            if ( ImguiProcessedEvent( &event ) ) {
                continue;
            }

            // saw any input at all, begin a short timer to allow inputs to accumulate
            input_latency = SDL_GetTicks();

            switch ( event.type ) {
            case SDL_MOUSEBUTTONDOWN:
                HandleMouseClick( event.button.x, event.button.y, &event.button );
                break;

            case SDL_MOUSEMOTION:
                HandleMouseMove( event.motion.x, event.motion.y );
                break;

            case SDL_MOUSEWHEEL:
                HandleMouseWheel( event.wheel.x, event.wheel.y );
                break;

            case SDL_QUIT:
                loop = false;
                break;

            case SDL_KEYDOWN:
                loop = HandleKeyPress( &event.key.keysym );
                break;

            default:
                break;
            }
        }
    }


    void my_stbtt_initfont(void)
    {
        // ALLOCATED
        unsigned char * ttf_buffer = new unsigned char[ 8 << 20 ];
        unsigned char * temp_bitmap = new unsigned char[ 512 * 512 ];
        cdata = (stbtt_bakedchar *) calloc( sizeof(stbtt_bakedchar), 96 );
        ftex = new GLuint[1];

        FILE * fh = fopen( FONT_FILENAME, "rb" );
        if ( !fh ) {
            printf( "couldn't open font: \"%s\"\n", FONT_FILENAME );
            return;
        }
        size_t whatever = fread( ttf_buffer, 1, 1<<20, fh );
        (void)whatever;
        stbtt_BakeFontBitmap(ttf_buffer, 0, font_px, temp_bitmap, font_scale, font_scale, 32, 96, cdata); // no guarantee this fits!

        // can free ttf_buffer at this point
        fclose( fh );
        delete[] ttf_buffer;

        glGenTextures(1, ftex);
        glBindTexture(GL_TEXTURE_2D, *ftex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, font_scale, font_scale, 0, GL_ALPHA, GL_UNSIGNED_BYTE, temp_bitmap);
        //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

        // can free temp_bitmap at this point
        delete[] temp_bitmap;
    }

    void my_stbtt_cleanup() {
        if ( cdata )
            free( cdata );
        if ( ftex )
            delete[] ftex;
    }

    void my_stbtt_print(float x, float y, char *text )
    {
        glColor4f( 1.0f, 1.0f, 1.0f, 1.0f );

        // assume orthographic projection with units = screen pixels, origin at top left
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, *ftex);
        glBegin(GL_QUADS);
        while (*text) {
            if (*text >= 32 && *text < 128) {
                stbtt_aligned_quad q;
                stbtt_GetBakedQuad(cdata, font_scale, font_scale, *text-32, &x, &y, &q, 1);//1=opengl & d3d10+,0=d3d9

                const stbtt_bakedchar *b = cdata + *text-32;
                int y_off = b->yoff;
                if ( *text == '-' || *text == '=' ) y_off = y_off + y_off/2;

                glTexCoord2f(q.s0,q.t1); glVertex2f(q.x0,q.y0-y_off);
                glTexCoord2f(q.s0,q.t0); glVertex2f(q.x0,q.y1-y_off);
                glTexCoord2f(q.s1,q.t0); glVertex2f(q.x1,q.y1-y_off);
                glTexCoord2f(q.s1,q.t1); glVertex2f(q.x1,q.y0-y_off);
            }
            ++text;
        }
        glEnd();
        glDisable(GL_TEXTURE_2D);
    }

    void DrawMessages() {
        char buf[256];
        double * lookAt = mandelbrot.GetLookAt();
        double * zoom = mandelbrot.GetZoom();
        sprintf( buf, "loc = %.3lf, %.3lf, zoom: %.5le, %.5le", lookAt[0], lookAt[1], zoom[0], zoom[1] );
        my_stbtt_print( windowWidth - 400, windowHeight - font_px, buf );

        if ( message_p ) {
            my_stbtt_print( windowWidth/2-70, windowHeight - font_px - 240, message_p );
        }

        double * hsv = mandelbrot.GetHSV();
        sprintf( buf, "H: %.2f, S: %.2f, V: %.2f", hsv[0], hsv[1], hsv[2] );
        my_stbtt_print( windowWidth - 200, 2, buf );

        if ( fileMessage.IsAlive() ) {
            my_stbtt_print( fileMessage.x, fileMessage.y, fileMessage.c_str() );
        }

        sprintf( buf, "%u.%u sec to compute", mandelbrot.GetComputeTime() / 1000, ( mandelbrot.GetComputeTime() % 1000 ) / 100 );
        my_stbtt_print( 2, 2, buf );
    }

    void ImguiDrawOverlay() {
        // Start the ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame( mainWindow );
        ImGui::NewFrame();


        // control panel window
        if ( gui_show_control_panel ) {
            ImGui::Begin( "Control Panel", &gui_show_control_panel );
            ImGui::Text( "Mandelbrot variables and configuration\n" );

            //-----------------
            ImGui::Text( " " );
            ImGui::Separator();
            ImGui::Text( " " );

            // static stats
            {
                int lcol = 190;

                ImGui::Columns( 2, "res", false );
                ImGui::SetColumnWidth( 0, lcol );
                ImGui::Text( "Window Resolution:" );
                ImGui::NextColumn();
                ImGui::Text( "%d, %d", windowWidth, windowHeight );
                ImGui::NextColumn();

                ImGui::Columns( 2, "mousecoord", false );
                ImGui::SetColumnWidth( 0, lcol );
                ImGui::Text( "Mouse Coordinates:" );
                ImGui::NextColumn();
                ImGui::Text( "%d, %d", currentMouseLocation[0], currentMouseLocation[1] );
                ImGui::NextColumn();

                ImGui::Columns( 2, "time2compute", false );
                ImGui::SetColumnWidth( 0, lcol );
                ImGui::Text( "Frame compute time:" );
                ImGui::NextColumn();
                ImGui::Text( "%u.%u sec", mandelbrot.GetComputeTime() / 1000, ( mandelbrot.GetComputeTime() % 1000 ) / 100 );
                ImGui::NextColumn();

                ImGui::Columns(1);
            }

            ImGui::Text( " " );

            // zoom and location
            {
                ImGui::Text( "Fundamental Characteristics:" );
                if ( ImGui::IsItemHovered() )
                    ImGui::SetTooltip("These 4 numbers are all you need to get back to this location.");

                double * la = mandelbrot.GetLookAt();
                double * zoom = mandelbrot.GetZoom();

                ImGui::Columns( 3, "center", false );
                ImGui::SetColumnWidth( 0, ImGui::GetWindowWidth() * 0.21f );
                ImGui::Text( "CENTER: " );
                if ( ImGui::IsItemHovered() )
                    ImGui::SetTooltip( "The center of the window location over the fractal" );
                ImGui::NextColumn();
                ImGui::SetColumnWidth( 1, ImGui::GetWindowWidth() * 0.36f );
                ImGui::InputScalar( "c.x",  ImGuiDataType_Double, &la[0], NULL );
                ImGui::NextColumn();
                ImGui::SetColumnWidth( 2, ImGui::GetWindowWidth() * 0.36f );
                ImGui::InputScalar( "c.y",  ImGuiDataType_Double, &la[1], NULL );
                ImGui::NextColumn();

                ImGui::Columns( 3, "zoom", false );
                ImGui::SetColumnWidth( 0, ImGui::GetWindowWidth() * 0.21f );
                ImGui::Text( "ZOOM: " );
                if ( ImGui::IsItemHovered() ) {
                    ImGui::BeginTooltip();
                    ImGui::Text( "The degree of magnification in the X and Y axises." );
                    ImGui::Text( "Another way to think about Zoom is that the number signifies" );
                    ImGui::Text( "half the measurement from one side of the viewport to the other" );
                    ImGui::Text( "of a single axis. Smaller viewport equals larger magnification." );
                    ImGui::EndTooltip();
                }
                ImGui::NextColumn();
                ImGui::SetColumnWidth( 1, ImGui::GetWindowWidth() * 0.36f );
                ImGui::InputScalar( "z.x",  ImGuiDataType_Double, &zoom[0], NULL );
                ImGui::NextColumn();
                ImGui::SetColumnWidth( 2, ImGui::GetWindowWidth() * 0.36f );
                ImGui::InputScalar( "z.y",  ImGuiDataType_Double, &zoom[1], NULL );
                ImGui::NextColumn();

                ImGui::Columns( 1 );
            }

            ImGui::Text( " " );

            // Color mode Combo dropdown
            {
                const char* items[] = { "1. RGB - red green blue", "2. RBG - red blue green", "3. GRB - green red blue", "4. GBR - green blue red", "5. BRG - blue red green", "6. BGR - blue green red" };
                MandelbrotV2::filter_t cFilter = mandelbrot.GetColorFilter();
                int item_current = (int) cFilter;
                ImGui::Combo( "Color Mode", &item_current, items, IM_ARRAYSIZE(items) );
                ImGui::SameLine(); ShowHelpMarker( "Pixels are calculated through various kinds of filters. Each filter places different emphasis on color groups which achieves vastly different effects." );
                if ( item_current != (int) cFilter ) {
                    mandelbrot.SetColorFilter( (MandelbrotV2::filter_t) item_current );
                }
            }

            ImGui::Text( " " );

            // HSV sliders
            {
                double * hue       ;
                double * saturation;
                double * vibrance  ;

                mandelbrot.GetHue( &hue );
                mandelbrot.GetSaturation( &saturation );
                mandelbrot.GetVibrance( &vibrance );

                double hlow = 10.0, hhigh = 1500.0;
                double slow = 0.1 , shigh = 1.0;
                double vlow = 1.0 , vhigh = 50.0;

                ImGui::SliderScalar( "Hue", ImGuiDataType_Double, hue, &hlow, &hhigh, "%.3lf", 1.0f );
                ImGui::SliderScalar( "Saturation", ImGuiDataType_Double, saturation, &slow, &shigh, "%.3lf", 1.0f );
                ImGui::SliderScalar( "Vibrance", ImGuiDataType_Double, vibrance, &vlow, &vhigh, "%.3lf", 1.0f );
            }

            ImGui::Text( " " );

            // controls to open other windows
            {
                ImGui::Checkbox( "Help Menu", &gui_show_help_window );
                ImGui::SameLine();
                ImGui::Checkbox( "ImGui Demo", &gui_show_demo_window );
                ImGui::SameLine();
                ImGui::Checkbox( "OpenGL Filtering", &gl_filter_on );
            }

            //-----------------
            ImGui::Text( " " );
            ImGui::Separator();
            ImGui::Text( " " );

            // buttons
            {
                if ( ImGui::Button( "Re-Center" ) )
                    gui_reset_location_requested = true;
                ImGui::SameLine();
                if ( ImGui::Button( "Reset All" ) )
                    gui_reset_all_requested = true;

                if ( ImGui::Button( "Screenshot" ) )
                    gui_screenshot_requested = true;
                ImGui::SameLine();
                if ( ImGui::Button( "\n   Draw   \n   " ) )
                    gui_redraw_requested = true;

                if ( ImGui::Button( "Quit" ) )
                    gui_quit_requested = true;
            }

            ImGui::End();
        } // control panel


        // WINDOW 2 - help window
        if ( gui_show_help_window )
        {
            struct cmd_t {
                const char * cmd;
                const char * desc;
            } cmds[] = {
                "ACTION or KEY", "DESCRIPTION",
                "", "",
                "Left-click", "moves center of view",
                "Mouse-wheel", "Zooms in or out",
                "G", "show/hide Control Panel Gui",
                "? or /", "show/hide this help window",
                "F", "toggle fullscreen",
                "H", "increase hue",
                "shift+H", "decrease hue",
                "S", "increase saturation",
                "shift+S", "decrease saturation",
                "V", "increase vibrance",
                "shift+V", "decrease vibrance",
                "M", "cycle color modes",
                "shift+M", "cycle color modes backwards",
                "O", "saves screenshot of the screen",
                "Q or ESCAPE", "quits/exits the program",
                "C", "changes or disables crosshair",
                "W", "Toggle OpenGl Filtering",
                nullptr, nullptr
            };

            ImGui::Begin( "Help", &gui_show_help_window );
            ImGui::Text( "Key Bindings" );
            ImGui::Text( " " );
            ImGui::Separator();
            cmd_t * hp = &cmds[0];

            int lcolwidth = 137;

            do {
                cmd_t & c = *hp++;
                if ( c.cmd ) {
                    if ( *c.cmd ) {
                        ImGui::Columns( 2, "key description", false );
                        ImGui::SetColumnWidth( 0, lcolwidth );
                        ImGui::Text( "(%s)", c.cmd );
                        ImGui::NextColumn();
                        ImGui::Text( "%s", c.desc );
                        ImGui::NextColumn();
                    } else {
                        ImGui::Columns( 1 );
                        ImGui::Separator();
                        ImGui::NextColumn();
                    }
                } else {
                    break;
                }
            } while ( 1 );

            ImGui::Separator();
            if (ImGui::Button("Close"))
                gui_show_help_window = false;
            ImGui::End();
        }

        // Window 3 - Show the ImGui demo window.
        if ( gui_show_demo_window )
        {
            ImGui::SetNextWindowPos( ImVec2(650, 20), ImGuiCond_FirstUseEver );
            ImGui::ShowDemoWindow( &gui_show_demo_window );
        }

        // DRAW
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData( ImGui::GetDrawData() );
    }

    void Render() {
        // init for drawing
        glMatrixMode( GL_MODELVIEW );
        glClear( GL_COLOR_BUFFER_BIT );
        SDL_GL_MakeCurrent( mainWindow, mainGLContext );

        // draw the computed mandelbrot
        DrawScreenTexture( textureHandle, screenTextureBuffer );

        // mouse line, crosshair and text goes behind gui
        //DrawMouseline();
        DrawCrosshairs();
        //DrawMessages();

        // draw the GUI
        ImguiDrawOverlay();

        // done drawing GL
        glFlush();  // flush all GL commands
        glFinish(); // block until all GL execution is complete

        SDL_GL_SwapWindow( mainWindow );
    }

    void Draw( bool recompute =false )
    {
        static char buf[32] = "RECOMPUTING. . . .";

        if ( recompute ) {
            message_p = buf;
            Render();

            mandelbrot.Compute( screenTextureBuffer );
        }

        message_p = nullptr;
        Render();
    }

    void Run() {
        // Clear buffer with black background. Same as: SDL_RenderClear(&renderer);
        glClear(GL_COLOR_BUFFER_BIT);
        SDL_GL_SwapWindow( mainWindow );

        if ( fullscreen_flag ) {
            SDL_SetWindowFullscreen( mainWindow, fullscreen_flag );
        }

        const int wait_for_input = 466; // num ms of pause to wait for until executing new instructions
        bool loop = true;
        int updateTime = 0;

        while ( loop && !gui_quit_requested )
        {
            PollEvents( loop, updateTime );

            // give the events a fraction of a second to accumulate before we block with a draw
            // this averts multiple draws for say a single flurry of mouse-wheel events or key-presses
            int t1 = SDL_GetTicks();
            if ( dataChanged && t1 - updateTime >= 0 && t1 - updateTime < wait_for_input ) {
                continue;
            }

            if ( gui_redraw_requested ) {
                dataChanged = true;
                gui_redraw_requested = false;
            }

            if ( gui_reset_all_requested ) {
                mandelbrot.ResetValues();
                ScalingFix();
                gui_reset_all_requested = false;
                dataChanged = true;
            }

            if ( gui_reset_location_requested ) {
                mandelbrot.ResetLocation();
                ScalingFix();
                gui_reset_location_requested = false;
                dataChanged = true;
            }

            Draw( dataChanged );
            dataChanged = false;

            if ( gui_screenshot_requested ) {
                gui_screenshot_requested = false;
                SaveScreenshot();
            }
        }
    }

    void Shutdown() {
        if ( !mainWindow || !screenTextureBuffer )
            return;

        // imgui
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplSDL2_Shutdown();
        ImGui::DestroyContext();

        //glDeleteShader(

        // Delete our OpengL context
        SDL_GL_DeleteContext( mainGLContext );

        // Destroy our window
        SDL_DestroyWindow( mainWindow );

        mainWindow = nullptr;

        // Shutdown SDL 2
        SDL_Quit();

        if ( screenTextureBuffer ) {
            free( screenTextureBuffer );
            screenTextureBuffer = nullptr;
        }
    }
}; // class Application


static int ParseArguments( int argc, char ** argv, const char * exename, int * w, int * h, bool & fullscreen ) {
    for ( int i = 1; i < argc; ++i ) {
        // set custom resoution
        if ( strcmp( argv[i], "-h" ) == 0 || strcmp( argv[i], "--help" ) == 0 ) {
            usage( exename );
            return -1;
        } else if ( strcmp( argv[i], "-f" ) == 0 ) {
            fullscreen = true;
            continue;
        } else if ( strcmp( argv[i], "-res" ) == 0 || strcmp( argv[i], "-xy" ) == 0 ) {
            if ( i+1 < argc ) {
                char buf[256];
                strncpy( buf, argv[i+1], sizeof(buf) );
				buf[sizeof(buf) - 1] = '\0';
                char * p = strchr( buf, 'x' );
                if ( !p ) {
                    usage( exename );
                    return -1;
                }
                *p++ = '\0';
                long int nx, ny;
                if ( !_getint( buf, &nx ) || !_getint( p, &ny ) ) {
                    usage( exename );
                    return -1;
                }
                *w = (signed int) nx;
                *h = (signed int) ny;
                ++i;
                continue;
            } else {
                usage( exename );
                return -1;
            }
        } else {
            usage( exename );
            return -1;
        }
    }
    return 0;
}

int main( int argc, char ** argv )
{
    int width = DEFAULT_WINDOW_WIDTH;
    int height = DEFAULT_WINDOW_HEIGHT;

    char * p = strrchr( argv[0], '/' );
    const char * exename = p ? p + 1 : argv[0];

    bool fullscreen = false;
    (void)fullscreen;

    // arguments
    int res = ParseArguments( argc, argv, exename, &width, &height, fullscreen );
    if ( res != 0 )
        return res;

    bool userHasSpecifiedResolution = ( width != DEFAULT_WINDOW_WIDTH ) ||
                                      ( height != DEFAULT_WINDOW_HEIGHT );

    // FIXME: should be singleton
    //  Application & app = Application::GetApplication();

    // declare our app
    Application app( width, height, fullscreen, !userHasSpecifiedResolution );

    // initialize
    if ( !app.Init() )
        return -1;

    // print usage hint
    if ( !userHasSpecifiedResolution ) {
        printf( "(hint: use '%s -xy <width>x<height>' to set custom resolution, -h for full help)\n", exename );
    }

    // display info
    app.Print_SDL_GL_Attributes();

    // run
    app.Run();

    return 0;
}
