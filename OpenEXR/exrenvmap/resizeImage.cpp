///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2007, Industrial Light & Magic, a division of Lucas
// Digital Ltd. LLC
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// *       Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// *       Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
// *       Neither the name of Industrial Light & Magic nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission. 
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------
//
//	resizeLatLong(), resizeCube() -- functions that resample
//	an environment map and convert it to latitude-longitude or
//	cube-face format.
//
//-----------------------------------------------------------------------------

#include <resizeImage.h>

#include "Iex.h"
#include <string.h>


using namespace std;
using namespace Imf;
using namespace Imath;


void
resizeLatLong (const EnvmapImage &image1,
	       EnvmapImage &image2,
	       const Box2i &image2DataWindow,
	       float filterRadius,
	       int numSamples)
{
    int w = image2DataWindow.max.x - image2DataWindow.min.x + 1;
    int h = image2DataWindow.max.y - image2DataWindow.min.y + 1;
    float radius = 0.5f * float(2 * M_PI) * filterRadius / w;

    image2.resize (ENVMAP_LATLONG, image2DataWindow);
    image2.clear ();

    Array2D<Rgba> &pixels = image2.pixels();

    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
	        V3f dir = LatLongMap::direction (image2DataWindow, V2f ((float) x, (float) y));
	        pixels[y][x] = image1.filteredLookup (dir, radius, numSamples);
        }
    }
}


namespace {
    
    V2f
    dirToPosLatLong (const Box2i &dataWindow, const V3f &dir)
    {
        return LatLongMap::pixelPosition (dataWindow, dir);
    }

    // values less than 1 cause ringing, but appear sharper.
    const float filterScale = 1.0f;

    inline float clean(float t)
    {
        const float EPSILON = .0000125f;
        if (fabs(t) < EPSILON)
            return 0.0f;
        return t;
    }
    
    inline float sinc(float x)
    {
        x = (x * float(M_PI));
        
        if ((x < 0.01f) && (x > -0.01f))
            return 1.0f + x * x * (-1.0f/6.0f + x * x * 1.0f/120.0f);
        
        return sinf(x) / x;
    }
    
    // see http://code.google.com/p/imageresampler/source/browse/trunk/resampler.cpp
    // for more filters
    
    const float lanczos3Support = 3.0f;
    float lanczosFilter(float t, float support)
    {
        if (t < 0.0f)
            t = -t;
        
        if (t < support)
            return clean(sinc(t) * sinc(t / support));
        else
            return (0.0f);
    }
    

    float calcFilterHalfWidth(float filterScale, float filterSupport, float inputSize, float outputSize)
    {
        // stretched half-width of filter
        float xscale = outputSize / inputSize;
        return (filterSupport / xscale) * filterScale;
    }
    
}



void
resizeCube (const EnvmapImage &srcImage,
            EnvmapImage &dstImage,
            const Box2i &dstImageDW,
            float filterRadius,
            int numSamples)
{
    dstImage.resize (ENVMAP_CUBE, dstImageDW);

    if (srcImage.type() == ENVMAP_CUBE && srcImage.dataWindow() == dstImageDW)
    {
        //
        // Special case - the input image is a cube-face environment
        // map with the same size as the output image.  We can copy
        // the input image without resampling.
        // 

        int w = dstImageDW.max.x - dstImageDW.min.x + 1;
        int h = dstImageDW.max.y - dstImageDW.min.y + 1;

        memcpy (&(dstImage.pixels()[0][0]),
                &(srcImage.pixels()[0][0]),
                sizeof (Rgba) * w * h);

        return;
    }
    
    const Box2i& srcImageDW = srcImage.dataWindow();
    int srcSof = CubeMap::sizeOfFace (srcImageDW);   // picks min of width and heigth/6
    int dstSof = CubeMap::sizeOfFace (dstImageDW);
    float radius = 1.5f * filterRadius / dstSof;
            
    Array2D<Rgba> &pixels = dstImage.pixels();

    if (srcImage.type() == ENVMAP_LATLONG)
    {
        for (int f = CUBEFACE_POS_X; f <= CUBEFACE_NEG_Z; ++f)
        {
            CubeMapFace face = CubeMapFace (f);
            
            for (int y = 0; y < dstSof; ++y)
            {
                for (int x = 0; x < dstSof; ++x)
                {
                    V3f dir = dstImage.direction(face, x, y);
                    V2f samplePos = dirToPosLatLong(srcImageDW, dir);
                    Imf::Rgba sample = srcImage.pixels()[int(floor(samplePos.y+0.5f))][int(floor(samplePos.x+0.5f))];
                    
                    dir = dstImage.direction(face, (x+1)%dstSof, y);
                    samplePos = V2f(dirToPosLatLong(srcImageDW, dir));
                    Imf::Rgba s2 = srcImage.pixels()[int(floor(samplePos.y+0.5f))][int(floor(samplePos.x+0.5f))];

                    dir = dstImage.direction(face, (x+1)%dstSof, (y+1)%dstSof);
                    samplePos = V2f(dirToPosLatLong(srcImageDW, dir));
                    Imf::Rgba s3 = srcImage.pixels()[int(floor(samplePos.y+0.5f))][int(floor(samplePos.x+0.5f))];
                    
                    dir = dstImage.direction(face, x, (y+1)%dstSof);
                    samplePos = V2f(dirToPosLatLong(srcImageDW, dir));
                    Imf::Rgba s4 = srcImage.pixels()[int(floor(samplePos.y+0.5f))][int(floor(samplePos.x+0.5f))];

                    sample.r += s2.r + s3.r + s4.r;
                    sample.g += s2.g + s3.g + s4.g;
                    sample.b += s2.b + s3.b + s4.b;
                    sample.a += s2.a + s3.a + s4.a;
                    sample.r *= 0.25f;
                    sample.g *= 0.25f;
                    sample.b *= 0.25f;
                    sample.a *= 0.25f;

                    V2f pos = dstImage.pixelPos(face, x, y);
                    pixels[int(floorf(pos.y+0.5f))][int(floorf(pos.x+0.5f))] = sample;
                }
            }
        }
        
        return;
    }
    
    if (srcImage.type() == ENVMAP_CUBE)
    {
        dstImage.precalcTables();
        for (int f = CUBEFACE_POS_X; f <= CUBEFACE_NEG_Z; ++f)
        {
            CubeMapFace face = CubeMapFace (f);
            
            for (int y = 0; y < dstSof; ++y)
            {
                for (int x = 0; x < dstSof; ++x)
                {
                    const bool convolve = true;
                    if (convolve)
                    {
                        V3f dstDir = dstImage.direction(f, x, y);
                        
                        V2f pos;
                        V3f sampleDir;
                        const float phongPower = 23;
                        float rTotal = 0;
                        float gTotal = 0;
                        float bTotal = 0;
                        float aTotal = 0;
                        double weightTotal = 0.0;
                        for (int f1 = CUBEFACE_POS_X; f1 <= CUBEFACE_NEG_Z; ++f1)
                        {
                            for (int y1 = 0; y1 < srcSof; ++y1)
                            {
                                for (int x1 = 0; x1 < srcSof; ++x1)
                                {
                                    sampleDir = srcImage.direction(f1, x1, y1);
                                    double weight = sampleDir ^ dstDir;
                                    if (weight <= 0)
                                        continue;
                                    
                                    weight = pow(weight, phongPower);
                                    weight *= srcImage.solidAngleWeight(f1, x1, y1);
                                    
                                    if (weight <= FLT_EPSILON)
                                        continue;
                                    
                                    weightTotal += weight;
                                    
                                    pos = srcImage.pixelPos(f1, x1, y1);
                                    Imf::Rgba sample = srcImage.pixels()[int(floor(pos.y+0.5f))][int(floor(pos.x+0.5f))];

                                    rTotal += sample.r * weight;
                                    gTotal += sample.g * weight;
                                    bTotal += sample.b * weight;
                                    aTotal += sample.a * weight;
                                }
                            }
                        }
                        double k = weightTotal > 0.0 ? 1.0 / weightTotal : 0.0;
                        Imf::Rgba sample(float(rTotal * k), float(gTotal * k), float(bTotal * k), float(aTotal * k));
                        pos = dstImage.pixelPos(f, x, y);
                        pixels[int(floorf(pos.y+0.5f))][int(floorf(pos.x+0.5f))] = sample;
                    }
                    else
                    {
                        CubeMapFace srcFace;
                        V2f pos;
                        V3f dir;

                        dir = dstImage.direction(face, x, y);
                        CubeMap::faceAndPixelPosition (dir, srcImageDW, srcFace, pos);
                        pos = srcImage.pixelPos(srcFace, floorf(pos.x+0.5f), floorf(pos.y+0.5f));
                        Imf::Rgba sample = srcImage.pixels()[int(floor(pos.y+0.5f))][int(floor(pos.x+0.5f))];
                        
                        dir = dstImage.direction(face, (x+1)%dstSof, y);
                        CubeMap::faceAndPixelPosition (dir, srcImageDW, srcFace, pos);
                        pos = srcImage.pixelPos(srcFace, floorf(pos.x+0.5f), floorf(pos.y+0.5f));
                        Imf::Rgba s2 = srcImage.pixels()[int(floor(pos.y+0.5f))][int(floor(pos.x+0.5f))];
                        
                        dir = dstImage.direction(face, (x+1)%dstSof, (y+1)%dstSof);
                        CubeMap::faceAndPixelPosition (dir, srcImageDW, srcFace, pos);
                        pos = srcImage.pixelPos(srcFace, floorf(pos.x+0.5f), floorf(pos.y+0.5f));
                        Imf::Rgba s3 = srcImage.pixels()[int(floor(pos.y+0.5f))][int(floor(pos.x+0.5f))];
                        
                        dir = dstImage.direction(face, x, (y+1)%dstSof);
                        CubeMap::faceAndPixelPosition (dir, srcImageDW, srcFace, pos);
                        pos = srcImage.pixelPos(srcFace, floorf(pos.x+0.5f), floorf(pos.y+0.5f));
                        Imf::Rgba s4 = srcImage.pixels()[int(floor(pos.y+0.5f))][int(floor(pos.x+0.5f))];
                        
                        sample.r += s2.r + s3.r + s4.r;
                        sample.g += s2.g + s3.g + s4.g;
                        sample.b += s2.b + s3.b + s4.b;
                        sample.a += s2.a + s3.a + s4.a;
                        sample.r *= 0.25f;
                        sample.g *= 0.25f;
                        sample.b *= 0.25f;
                        sample.a *= 0.25f;
                        
                        pos = dstImage.pixelPos(face, x, y);
                        pixels[int(floorf(pos.y+0.5f))][int(floorf(pos.x+0.5f))] = sample;
                    }
                }
            }
        }
        
        return;
    }

    //
    // Resample the input image for all other types of input maps
    //

    dstImage.clear ();

    for (int f = CUBEFACE_POS_X; f <= CUBEFACE_NEG_Z; ++f)
    {
        CubeMapFace face = CubeMapFace (f);

        for (int y = 0; y < dstSof; ++y)
        {
            for (int x = 0; x < dstSof; ++x)
            {
                V2f posInFace ((float) x, (float) y);

                V3f dir =
                    CubeMap::direction (face, dstImageDW, posInFace);

                V2f pos =
                    CubeMap::pixelPosition (face, dstImageDW, posInFace);

                pixels[int (pos.y + 0.5f)][int (pos.x + 0.5f)] =
                    srcImage.filteredLookup (dir, radius, numSamples);
            }
        }
    }
}
