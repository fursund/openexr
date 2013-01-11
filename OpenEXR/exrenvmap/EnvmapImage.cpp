///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2004, Industrial Light & Magic, a division of Lucas
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
//	class EnvmapImage
//
//-----------------------------------------------------------------------------

#include <EnvmapImage.h>
#include "ImathFun.h"


using namespace Imf;
using namespace Imath;


EnvmapImage::EnvmapImage():
    _type (ENVMAP_LATLONG),
    _dataWindow (V2i (0, 0), V2i (0, 0)),
    _pixels (1, 1),
    _solidAngleWeight(0),
    _directions(0),
    _positionsInFace(0),
    _sof(0)
{
    clear();
}

EnvmapImage::~EnvmapImage()
{
    delete[] _solidAngleWeight;
    delete[] _directions;
    delete[] _positionsInFace;
}

EnvmapImage::EnvmapImage (Envmap type, const Box2i &dataWindow):
    _type (type),
    _dataWindow (dataWindow),
    _pixels (dataWindow.max.y - dataWindow.min.y + 1,
             dataWindow.max.x - dataWindow.min.x + 1),
    _solidAngleWeight(0),
    _directions(0),
    _positionsInFace(0),
    _sof(0)
{
    clear();
}


void
EnvmapImage::resize (Envmap type, const Box2i &dataWindow)
{
    _pixels.resizeEraseUnsafe (dataWindow.max.y - dataWindow.min.y + 1,
                               dataWindow.max.x - dataWindow.min.x + 1);
    _type = type;
    _dataWindow = dataWindow;

    precalcTables();
    clear();
}


void
EnvmapImage::clear ()
{
    _sof = CubeMap::sizeOfFace (_dataWindow);   // picks min of width and heigth/6

    int w = _dataWindow.max.x - _dataWindow.min.x + 1;
    int h = _dataWindow.max.y - _dataWindow.min.y + 1;

    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            Rgba &p = _pixels[y][x];
            
            p.r = 0;
            p.g = 0;
            p.b = 0;
            p.a = 0;
        }
    }
}






namespace {

V2f
dirToPosLatLong (const Box2i &dataWindow, const V3f &dir)
{
    return LatLongMap::pixelPosition (dataWindow, dir);
}


V2f
dirToPosCube (const Box2i &dataWindow, const V3f &dir)
{
    CubeMapFace face;
    V2f posInFace;
    CubeMap::faceAndPixelPosition (dir, dataWindow, face, posInFace);
    return CubeMap::pixelPosition (face, dataWindow, posInFace);
}
    

} // namespace


Rgba
EnvmapImage::filteredLookup (V3f d, float r, int n) const
{
    //
    // Filtered environment map lookup: Take n by n point samples
    // from the environment map, clustered around direction d, and
    // combine the samples with a tent filter.
    //
    
    //
    // Depending on the type of map, pick an appropriate function
    // to convert 3D directions to 2D pixel poitions.
    //

    V2f (* dirToPos) (const Box2i &, const V3f &);

    if (_type == ENVMAP_LATLONG)
        dirToPos = dirToPosLatLong;
    else
        dirToPos = dirToPosCube;

    //
    // Pick two vectors, dx and dy, of length r, that are orthogonal
    // to the lookup direction, d, and to each other.
    //

    d.normalize();
    V3f dx, dy;

    if (abs (d.x) > 0.707f)
        dx = (d % V3f (0, 1, 0)).normalized() * r;
    else
        dx = (d % V3f (1, 0, 0)).normalized() * r;

    dy = (d % dx).normalized() * r;

    //
    // Take n by n point samples from the map, and add them up.
    // The directions for the point samples are all within the pyramid
    // defined by the vectors d-dy-dx, d-dy+dx, d+dy-dx, d+dy+dx.
    //

    float wt = 0;

    float cr = 0;
    float cg = 0;
    float cb = 0;
    float ca = 0;

    for (int y = 0; y < n; ++y)
    {
        float ry = float (2 * y + 2) / float (n + 1) - 1;
        float wy = 1 - abs (ry);
        V3f ddy (ry * dy);
        
        for (int x = 0; x < n; ++x)
        {
            float rx = float (2 * x + 2) / float (n + 1) - 1;
            float wx = 1 - abs (rx);
            V3f ddx (rx * dx);
            
            Rgba s = sample (dirToPos (_dataWindow, d + ddx + ddy));
            
            float w = wx * wy;
            wt += w;
            
            cr += s.r * w;
            cg += s.g * w;
            cb += s.b * w;
            ca += s.a * w;
        }
    }

    wt = 1 / wt;

    Rgba c;

    c.r = cr * wt;
    c.g = cg * wt;
    c.b = cb * wt;
    c.a = ca * wt;

    return c;
}


Rgba
EnvmapImage::sample (const V2f &pos) const
{
    //
    // Point-sample the environment map image at 2D position pos.
    // Interpolate bilinearly between the four nearest pixels.
    //

    int x1 = Imath::floor (pos.x);
    int x2 = x1 + 1;
    float sx = x2 - pos.x;
    float tx = 1 - sx;

    x1 = clamp (x1, _dataWindow.min.x, _dataWindow.max.x) - _dataWindow.min.x;
    x2 = clamp (x2, _dataWindow.min.x, _dataWindow.max.x) - _dataWindow.min.x;

    int y1 = Imath::floor (pos.y);
    int y2 = y1 + 1;
    float sy = y2 - pos.y;
    float ty = 1 - sy;

    y1 = clamp (y1, _dataWindow.min.y, _dataWindow.max.y) - _dataWindow.min.y;
    y2 = clamp (y2, _dataWindow.min.y, _dataWindow.max.y) - _dataWindow.min.y;

    Rgba p11 = _pixels[y1][x1];
    Rgba p12 = _pixels[y1][x2];
    Rgba p21 = _pixels[y2][x1];
    Rgba p22 = _pixels[y2][x2];

    Rgba p;
    p.r = (p11.r * sx + p12.r * tx) * sy + (p21.r * sx + p22.r * tx) * ty;
    p.g = (p11.g * sx + p12.g * tx) * sy + (p21.g * sx + p22.g * tx) * ty;
    p.b = (p11.b * sx + p12.b * tx) * sy + (p21.b * sx + p22.b * tx) * ty;
    p.a = (p11.a * sx + p12.a * tx) * sy + (p21.a * sx + p22.a * tx) * ty;

    return p;
}


namespace {

    
    // faces are ordered 0-5, as XYZ, posneg
    V3f cardinalDirPerFace(int face, int& ix, int& iy, int& iz)
    {
        V3f faceDir;
        switch (face)
        {
            case CUBEFACE_POS_X:
                faceDir = V3f (1, 0, 0);
                ix = 0;
                iy = 1;
                iz = 2;
                break;
                
            case CUBEFACE_NEG_X:
                faceDir = V3f (-1, 0, 0);
                ix = 0;
                iy = 1;
                iz = 2;
                break;
                
            case CUBEFACE_POS_Y:
                faceDir = V3f (0, 1, 0);
                ix = 1;
                iy = 0;
                iz = 2;
                break;
                
            case CUBEFACE_NEG_Y:
                faceDir = V3f (0, -1, 0);
                ix = 1;
                iy = 0;
                iz = 2;
                break;
                
            case CUBEFACE_POS_Z:
                faceDir = V3f (0, 0, 1);
                ix = 2;
                iy = 0;
                iz = 1;
                break;
                
            case CUBEFACE_NEG_Z:
                faceDir = V3f (0, 0, -1);
                ix = 2;
                iy = 0;
                iz = 1;
                break;
        }
        return faceDir;
    }
    
    inline double
    sqr (double x)
    {
        return x * x;
    }


} // anon


void
EnvmapImage::precalcTables() const
{
    _sof = CubeMap::sizeOfFace (_dataWindow);
    _solidAngleWeight = new double[6 * _sof * _sof];
    _directions = new Imath::V3f[6 * _sof * _sof];
    _positionsInFace = new Imath::V2f[6 * _sof * _sof];
    
    for (int f = CUBEFACE_POS_X; f <= CUBEFACE_NEG_Z; ++f)
    {
        int ix = 0, iy = 0, iz = 0;
        CubeMapFace face = CubeMapFace (f);
        V3f faceDir = cardinalDirPerFace(f, ix, iy, iz);
        
        for (int y = 0; y < _sof; ++y)
        {
            bool yEdge = (y == 0 || y == _sof - 1);
            
            for (int x = 0; x < _sof; ++x)
            {
                bool xEdge = (x == 0 || x == _sof - 1);
                
                int index = weightIndex(f, x, y);
                
                double weight;
                V2f posInFace ((float) x, (float) y);
                V2f pos = CubeMap::pixelPosition (face, _dataWindow, posInFace);
                _positionsInFace[index] = pos;
                
                V3f dir = CubeMap::direction (face, _dataWindow, posInFace).normalized();
                _directions[index] = dir;
                
                //
                // The solid angle subtended by pixel (x,y), as seen
                // from the center of the cube, is proportional to the
                // square of the distance of the pixel from the center
                // of the cube and proportional to the dot product of
                // the viewing direction and the normal of the cube
                // face that contains the pixel.
                //
                
                weight = (dir ^ faceDir) *
                         (sqr (dir[iy] / dir[ix]) + sqr (dir[iz] / dir[ix]) + 1);
                
                //
                // Pixels at the edges and corners of the
                // cube are duplicated; we must adjust the
                // pixel weights accordingly.
                //
                
                if (xEdge && yEdge)
                    weight /= 3;
                else if (xEdge || yEdge)
                    weight /= 2;
                
                _solidAngleWeight[index] = weight;
            }
        }
    }
}


