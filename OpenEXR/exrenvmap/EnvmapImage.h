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

#ifndef INCLUDED_ENVMAP_IMAGE_H
#define INCLUDED_ENVMAP_IMAGE_H

//-----------------------------------------------------------------------------
//
//	class EnvmapImage
//
//-----------------------------------------------------------------------------

#include <ImfArray.h>
#include <ImfRgba.h>
#include <ImfEnvmap.h>
#include "ImathBox.h"


class EnvmapImage
{
  public:

    EnvmapImage ();
    EnvmapImage (Imf::Envmap type, const Imath::Box2i &dataWindow);
    ~EnvmapImage ();
    
    void resize (Imf::Envmap type,
                 const Imath::Box2i &dataWindow);
    
    void clear ();
    
    Imf::Envmap         type () const { return _type; }
    const Imath::Box2i& dataWindow () const { return _dataWindow; }
    
    
    Imf::Array2D<Imf::Rgba> &       pixels () { return _pixels; }
    const Imf::Array2D<Imf::Rgba> & pixels () const { return _pixels; }
    
    Imf::Rgba sample (const Imath::V2f &pos) const;
    Imf::Rgba filteredLookup (Imath::V3f direction,
                              float radius,
                              int numSamples) const;

    void precalcTables() const;
    
    inline int weightIndex(int f, int x, int y) const
    {
        return (f * _sof * _sof) + (y * _sof) + x;
    }

    inline Imath::V3f direction(int face, int x, int y) const
    {
        int index = weightIndex(face, x, y);
        return _directions[index];
    }

    inline double solidAngleWeight(int face, int x, int y) const
    {
        int index = weightIndex(face, x, y);
        return _solidAngleWeight[index];
    }
    
    inline Imath::V2f pixelPos(int face, int x, int y) const
    {
        int index = weightIndex(face, x, y);
        return _positionsInFace[index];
    }

    
  private:
    Imf::Envmap             _type;
    Imath::Box2i            _dataWindow;
    Imf::Array2D<Imf::Rgba> _pixels;
    
    mutable double*         _solidAngleWeight;
    mutable Imath::V3f*     _directions;
    mutable Imath::V2f*     _positionsInFace;
    mutable int             _sof;
};


#endif
