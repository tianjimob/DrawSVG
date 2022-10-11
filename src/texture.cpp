#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  for (size_t i = 1; i < tex.mipmap.size(); ++i) {
      MipLevel& mip = tex.mipmap[i];
      MipLevel& mip_prev = tex.mipmap[i - 1];
      auto& texels = mip_prev.texels;
      for (size_t y = 0; y < mip.height; ++y) {
          for (size_t x = 0; x < mip.width; ++x) {
              int index_prev = 2 * y * mip_prev.width + 2 * x;
              Color tl, tr, bl, br;
              uint8_to_float(&tl.r, &texels[4 * index_prev]);
              uint8_to_float(&tr.r, &texels[4 * (index_prev + 1)]);
              uint8_to_float(&bl.r, &texels[4 * (index_prev + mip_prev.width)]);
              uint8_to_float(&br.r, &texels[4 * (index_prev + mip_prev.width + 1)]);
              Color c = (tl + tr + bl + br) * 0.25;
              int index = y * mip.width + x;
              float_to_uint8(&mip.texels[4 * index], &c.r);
          }
      }
  }
}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 6: Implement nearest neighbour interpolation
  float u_p = u * tex.width - 0.5;
  float v_p = v * tex.height - 0.5;
  int u0 = std::floor(u_p);
  int v0 = std::floor(v_p);
  int u1 = u0 + 1;
  int v1 = v0 + 1;
  float u0_weight = u1 - u_p;
  float v0_weight = v1 - v_p;
  MipLevel& mipLevel = tex.mipmap[level];
  size_t width = mipLevel.width;
  if (u0_weight >= 0.5) {
      if (v0_weight >= 0.5)
          return Color{ mipLevel.texels[4 * (u0 + v0 * width)] / 255.f, mipLevel.texels[4 * (u0 + v0 * width) + 1] / 255.f, mipLevel.texels[4 * (u0 + v0 * width) + 2] / 255.f, mipLevel.texels[4 * (u0 + v0 * width) + 3] / 255.f };
      else
          return Color{ mipLevel.texels[4 * (u0 + v1 * width)] / 255.f, mipLevel.texels[4 * (u0 + v1 * width) + 1] / 255.f, mipLevel.texels[4 * (u0 + v1 * width) + 2] / 255.f, mipLevel.texels[4 * (u0 + v1 * width) + 3] / 255.f };
  }
  else {
      if(v0_weight >= 0.5)
          return Color{ mipLevel.texels[4 * (u1 + v0 * width)] / 255.f, mipLevel.texels[4 * (u1 + v0 * width) + 1] / 255.f, mipLevel.texels[4 * (u1 + v0 * width) + 2] / 255.f, mipLevel.texels[4 * (u1 + v0 * width) + 3] / 255.f };
      else
          return Color{ mipLevel.texels[4 * (u1 + v1 * width)] / 255.f, mipLevel.texels[4 * (u1 + v1 * width) + 1] / 255.f, mipLevel.texels[4 * (u1 + v1 * width) + 2] / 255.f, mipLevel.texels[4 * (u1 + v1 * width) + 3] / 255.f };
  }
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering
  MipLevel& mipLevel = tex.mipmap[level];
  float u_p = u * mipLevel.width - 0.5;
  float v_p = v * mipLevel.height - 0.5;
  int u0 = std::floor(u_p);
  int v0 = std::floor(v_p);
  int u1 = u0 + 1;
  int v1 = v0 + 1;
  float u0_weight = u1 - u_p;
  float v0_weight = v1 - v_p;
  size_t width = mipLevel.width;
  if (u0 >= 0 && u0 < width - 1 && v0 >= 0 && v0 < mipLevel.height - 1)
  {
      Color u0_v0 = { mipLevel.texels[4 * (u0 + v0 * width)] / 255.f, mipLevel.texels[4 * (u0 + v0 * width) + 1] / 255.f, mipLevel.texels[4 * (u0 + v0 * width) + 2] / 255.f, mipLevel.texels[4 * (u0 + v0 * width) + 3] / 255.f };
      Color u0_v1 = { mipLevel.texels[4 * (u0 + v1 * width)] / 255.f, mipLevel.texels[4 * (u0 + v1 * width) + 1] / 255.f, mipLevel.texels[4 * (u0 + v1 * width) + 2] / 255.f, mipLevel.texels[4 * (u0 + v1 * width) + 3] / 255.f };
      Color u1_v0 = { mipLevel.texels[4 * (u1 + v0 * width)] / 255.f, mipLevel.texels[4 * (u1 + v0 * width) + 1] / 255.f, mipLevel.texels[4 * (u1 + v0 * width) + 2] / 255.f, mipLevel.texels[4 * (u1 + v0 * width) + 3] / 255.f };
      Color u1_v1 = { mipLevel.texels[4 * (u1 + v1 * width)] / 255.f, mipLevel.texels[4 * (u1 + v1 * width) + 1] / 255.f, mipLevel.texels[4 * (u1 + v1 * width) + 2] / 255.f, mipLevel.texels[4 * (u1 + v1 * width) + 3] / 255.f };
      return u0_weight * v0_weight * u0_v0 + u0_weight * (1 - v0_weight) * u0_v1 + (1 - u0_weight) * v0_weight * u1_v0 + (1 - u0_weight) * (1 - v0_weight) * u1_v1;
  }
  else if (u0 >= 0 && u0 < width && v0 >= 0 && v0 < mipLevel.height)
      return { mipLevel.texels[4 * (u0 + v0 * width)] / 255.f, mipLevel.texels[4 * (u0 + v0 * width) + 1] / 255.f, mipLevel.texels[4 * (u0 + v0 * width) + 2] / 255.f, mipLevel.texels[4 * (u0 + v0 * width) + 3] / 255.f };
  else
      return Color::White;
}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering
  auto k0 = std::floor(u_scale);
  auto k1 = k0 + 1;
  auto w0 = k1 - u_scale;
  auto w1 = 1 - w0;
  k0 = k0 < 0 ? 0 : k0;
  k0 = k0 > tex.mipmap.size() - 1 ? tex.mipmap.size() - 1 : k0;
  k1 = k1 < 0 ? 0 : k1;
  k1 = k1 > tex.mipmap.size() - 1 ? tex.mipmap.size() - 1 : k1;
  Color c0 = sample_bilinear(tex, u, v, k0);
  Color c1 = sample_bilinear(tex, u, v,k1);
  return w0 * c0 + w1 * c1;
}

} // namespace CMU462
