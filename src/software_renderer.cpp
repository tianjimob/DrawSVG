#include "software_renderer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "triangulation.h"

using namespace std;

namespace CMU462 {

// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg(SVG& svg)
{
#if defined(SSAA) || defined(MLAA)
    clear_my_target();
#endif
    // set top level transformation
    transformation = svg_2_screen;

    // draw all elements
    for (size_t i = 0; i < svg.elements.size(); ++i) {
        draw_element(svg.elements[i]);
    }

    // draw canvas outline
    Vector2D a = transform(Vector2D(0, 0));
    a.x--;
    a.y--;
    Vector2D b = transform(Vector2D(svg.width, 0));
    b.x++;
    b.y--;
    Vector2D c = transform(Vector2D(0, svg.height));
    c.x--;
    c.y++;
    Vector2D d = transform(Vector2D(svg.width, svg.height));
    d.x++;
    d.y++;

    rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
    rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
    rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
    rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

    // resolve and send to render target
    resolve();
}

void SoftwareRendererImp::set_sample_rate(size_t sample_rate)
{

    // Task 4:
    // You may want to modify this for supersampling support
#ifdef SSAA
    supersample_target.resize(4 * sample_rate * sample_rate * target_w * target_h);
#elif defined(MLAA)
    mlaa_sample_target.resize(4 * target_w * target_h);
    mlaa = !mlaa;
#endif
    this->sample_rate = sample_rate;
}

void SoftwareRendererImp::set_render_target(unsigned char* render_target,
    size_t width, size_t height)
{

    // Task 4:
    // You may want to modify this for supersampling support
#ifdef SSAA
    supersample_target.resize(4 * sample_rate * sample_rate * width * height);
#elif defined(MLAA)
    mlaa_sample_target.resize(4 * width * width);
#endif
#if defined(SSAA) || defined(MLAA)
    render_target_tmp.resize(4 * width * height);
#endif
    this->render_target = render_target;
    this->target_w = width;
    this->target_h = height;
}

void SoftwareRendererImp::draw_element(SVGElement* element)
{

    // Task 5 (part 1):
    // Modify this to implement the transformation stack

    switch (element->type) {
    case POINT:
        draw_point(static_cast<Point&>(*element));
        break;
    case LINE:
        draw_line(static_cast<Line&>(*element));
        break;
    case POLYLINE:
        draw_polyline(static_cast<Polyline&>(*element));
        break;
    case RECT:
        draw_rect(static_cast<Rect&>(*element));
        break;
    case POLYGON:
        draw_polygon(static_cast<Polygon&>(*element));
        break;
    case ELLIPSE:
        draw_ellipse(static_cast<Ellipse&>(*element));
        break;
    case IMAGE:
        draw_image(static_cast<Image&>(*element));
        break;
    case GROUP:
        draw_group(static_cast<Group&>(*element));
        break;
    default:
        break;
    }
}

// Primitive Drawing //

void SoftwareRendererImp::draw_point(Point& point)
{

    Vector2D p = transform(point.position);
    rasterize_point(p.x, p.y, point.style.fillColor);
}

void SoftwareRendererImp::draw_line(Line& line)
{

    Vector2D p0 = transform(line.from);
    Vector2D p1 = transform(line.to);
    rasterize_line(p0.x, p0.y, p1.x, p1.y, line.style.strokeColor);
}

void SoftwareRendererImp::draw_polyline(Polyline& polyline)
{

    Color c = polyline.style.strokeColor;

    if (c.a != 0) {
        int nPoints = polyline.points.size();
        for (int i = 0; i < nPoints - 1; i++) {
            Vector2D p0 = transform(polyline.points[(i + 0) % nPoints]);
            Vector2D p1 = transform(polyline.points[(i + 1) % nPoints]);
            rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
        }
    }
}

void SoftwareRendererImp::draw_rect(Rect& rect)
{

    Color c;

    // draw as two triangles
    float x = rect.position.x;
    float y = rect.position.y;
    float w = rect.dimension.x;
    float h = rect.dimension.y;

    Vector2D p0 = transform(Vector2D(x, y));
    Vector2D p1 = transform(Vector2D(x + w, y));
    Vector2D p2 = transform(Vector2D(x, y + h));
    Vector2D p3 = transform(Vector2D(x + w, y + h));

    // draw fill
    c = rect.style.fillColor;
    if (c.a != 0) {
        rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
        rasterize_triangle(p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c);
    }

    // draw outline
    c = rect.style.strokeColor;
    if (c.a != 0) {
        rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
        rasterize_line(p1.x, p1.y, p3.x, p3.y, c);
        rasterize_line(p3.x, p3.y, p2.x, p2.y, c);
        rasterize_line(p2.x, p2.y, p0.x, p0.y, c);
    }
}

void SoftwareRendererImp::draw_polygon(Polygon& polygon)
{

    Color c;

    // draw fill
    c = polygon.style.fillColor;
    if (c.a != 0) {

        // triangulate
        vector<Vector2D> triangles;
        triangulate(polygon, triangles);

        // draw as triangles
        for (size_t i = 0; i < triangles.size(); i += 3) {
            Vector2D p0 = transform(triangles[i + 0]);
            Vector2D p1 = transform(triangles[i + 1]);
            Vector2D p2 = transform(triangles[i + 2]);
            rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
        }
    }

    // draw outline
    c = polygon.style.strokeColor;
    if (c.a != 0) {
        int nPoints = polygon.points.size();
        for (int i = 0; i < nPoints; i++) {
            Vector2D p0 = transform(polygon.points[(i + 0) % nPoints]);
            Vector2D p1 = transform(polygon.points[(i + 1) % nPoints]);
            rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
        }
    }
}

void SoftwareRendererImp::draw_ellipse(Ellipse& ellipse)
{

    // Extra credit
}

void SoftwareRendererImp::draw_image(Image& image)
{

    Vector2D p0 = transform(image.position);
    Vector2D p1 = transform(image.position + image.dimension);

    rasterize_image(p0.x, p0.y, p1.x, p1.y, image.tex);
}

void SoftwareRendererImp::draw_group(Group& group)
{

    for (size_t i = 0; i < group.elements.size(); ++i) {
        draw_element(group.elements[i]);
    }
}

// Rasterization //

inline static void render_point(int x, int y, Color& color,
    unsigned char* render_target, size_t target_w)
{
    render_target[4 * (x + y * target_w)] = (uint8_t)(color.r * 255);
    render_target[4 * (x + y * target_w) + 1] = (uint8_t)(color.g * 255);
    render_target[4 * (x + y * target_w) + 2] = (uint8_t)(color.b * 255);
    render_target[4 * (x + y * target_w) + 3] = (uint8_t)(color.a * 255);
}

// The input arguments in the rasterization functions
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point(float x, float y, Color color)
{

    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= target_w)
        return;
    if (sy < 0 || sy >= target_h)
        return;

        // fill sample - NOT doing alpha blending!
#if defined(SSAA) || defined(MLAA)
    render_point(x, y, color, &render_target_tmp[0], target_w);
#else
    render_point(x, y, color, render_target, target_w);
#endif
}

inline static float implicitLine(float x0, float y0, float x1, float y1,
    float x, float y)
{
    return (y0 - y1) * x + (x1 - x0) * y + x0 * y1 - x1 * y0;
}

void SoftwareRendererImp::rasterize_line(float x0, float y0, float x1, float y1,
    Color color)
{

    // Task 2:
    // Implement line rasterization
    if (x0 > x1) {
        swap(x0, x1);
        swap(y0, y1);
    }

#ifdef LINE_DRAWING_MIDPOINT
    // midpoint algorithm
    float slope = (y1 - y0) / (x1 - x0);
    if (slope > 0.f && slope <= 1.f) {
        float d = implicitLine(x0, y0, x1, y1, x0 + 1.f, y0 + 0.5f);
        auto deltaHit = x1 - x0 + y0 - y1;
        auto deltaMiss = y0 - y1;

        int y = (int)floor(y0);
        x0 = (int)floor(x0);
        x1 = (int)floor(x1);

        for (int x = x0; x <= x1; ++x) {
#if defined(SSAA) || defined(MLAA)
            render_point(x, y, color, &render_target_tmp[0], target_w);
#else
            render_point(x, y, color, render_target, target_w);
#endif
            if (d < 0.f) {
                y += 1;
                d += deltaHit;
            } else {
                d += deltaMiss;
            }
        }
    } else if (slope > 1.f) {
        float d = implicitLine(x0, y0, x1, y1, x0 + 0.5f, y0 + 1.f);
        auto deltaHit = x1 - x0 + y0 - y1;
        auto deltaMiss = x1 - x0;

        int x = (int)floor(x0);
        y0 = (int)floor(y0);
        y1 = (int)floor(y1);

        for (int y = y0; y <= y1; ++y) {
#if defined(SSAA) || defined(MLAA)
            render_point(x, y, color, &render_target_tmp[0], target_w);
#else
            render_point(x, y, color, render_target, target_w);
#endif
            if (d > 0.f) {
                x += 1;
                d += deltaHit;
            } else {
                d += deltaMiss;
            }
        }
    } else if (slope >= -1.f && slope <= 0.f) {
        float d = implicitLine(x0, y0, x1, y1, x0 + 1.f, y0 - 0.5f);
        auto deltaHit = x0 - x1 + y0 - y1;
        auto deltaMiss = y0 - y1;

        int y = (int)floor(y0);
        x0 = (int)floor(x0);
        x1 = (int)floor(x1);

        for (int x = x0; x <= x1; ++x) {
#if defined(SSAA) || defined(MLAA)
            render_point(x, y, color, &render_target_tmp[0], target_w);
#else
            render_point(x, y, color, render_target, target_w);
#endif
            if (d > 0.f) {
                y -= 1;
                d += deltaHit;
            } else {
                d += deltaMiss;
            }
        }
    } else if (slope < -1.f) {
        float d = implicitLine(x0, y0, x1, y1, x0 + 0.5f, y0 - 1.f);
        auto deltaHit = (float)(x0 - x1 + y0 - y1);
        auto deltaMiss = (float)(x0 - x1);

        int x = (int)floor(x0);
        y0 = (int)floor(y0);
        y1 = (int)floor(y1);

        for (int y = y0; y >= y1; --y) {
#if defined(SSAA) || defined(MLAA)
            render_point(x, y, color, &render_target_tmp[0], target_w);
#else
            render_point(x, y, color, render_target, target_w);
#endif
            if (d < 0.f) {
                x += 1;
                d += deltaHit;
            } else {
                d += deltaMiss;
            }
        }
    }

#elif defined(LINE_DRWAING_BRESENHAM)
    // bresenham algorithm
    int sx0 = (int)floor(x0);
    int sy0 = (int)floor(y0);
    int sx1 = (int)floor(x1);
    int sy1 = (int)floor(y1);
    float m = (y1 - y0) / (x1 - x0);
    int dx = x1 - x0;
    int dy = y1 - y0;
    if (m >= 0 && m <= 1) {
        int y = y0;
        int eps = 0;
        for (int x = sx0; x < sx1; ++x) {
#if defined(SSAA) || defined(MLAA)
            render_point(x, y, color, &render_target_tmp[0], target_w);
#else
            render_point(x, y, color, render_target, target_w);
#endif
            eps += dy;
            if ((eps << 1) > dx) {
                ++y;
                eps -= dx;
            }
        }
    } else if (m > 1) {
        int x = x0;
        int eps = 0;
        for (int y = sy0; y < sy1; ++y) {
#if defined(SSAA) || defined(MLAA)
            render_point(x, y, color, &render_target_tmp[0], target_w);
#else
            render_point(x, y, color, render_target, target_w);
#endif
            eps += dx;
            if ((eps << 1) > dy) {
                ++x;
                eps -= dy;
            }
        }
    } else if (m >= -1 && m < 0) {
        int y = y0;
        int eps = 0;
        for (int x = sx0; x < sx1; ++x) {
#if defined(SSAA) || defined(MLAA)
            render_point(x, y, color, &render_target_tmp[0], target_w);
#else
            render_point(x, y, color, render_target, target_w);
#endif
            eps += dy;
            if ((eps << 1) < -dx) {
                --y;
                eps += dx;
            }
        }
    } else {
        int x = x0;
        int eps = 0;
        for (int y = sy0; y > sy1; --y) {
#if defined(SSAA) || defined(MLAA)
            render_point(x, y, color, &render_target_tmp[0], target_w);
#else
            render_point(x, y, color, render_target, target_w);
#endif
            eps -= dx;
            if ((eps << 1) < dy) {
                ++x;
                eps -= dy;
            }
        }
    }
#endif
}

void SoftwareRendererImp::rasterize_triangle(float x0, float y0, float x1,
    float y1, float x2, float y2,
    Color color)
{
    // Task 3:
    // Implement triangle rasterization
#ifdef SSAA
    x0 *= sample_rate;
    y0 *= sample_rate;
    x1 *= sample_rate;
    y1 *= sample_rate;
    x2 *= sample_rate;
    y2 *= sample_rate;
#endif
    auto [xmin, xmax] = std::minmax({ x0, x1, x2 });
    auto [ymin, ymax] = std::minmax({ y0, y1, y2 });
    xmin = std::floor(xmin);
    xmax = std::ceil(xmax);
    ymin = std::floor(ymin);
    ymax = std::ceil(ymax);
    float dxAlpha = y1 - y2;
    float dyAlpha = x2 - x1;
    float dxBeta = y0 - y2;
    float dyBeta = x2 - x0;
    float dxGamma = y0 - y1;
    float dyGamma = x1 - x0;
    float numeratorAlpha = dxAlpha * (xmin + 0.5) + dyAlpha * (ymin + 0.5) + x1 * y2 - x2 * y1;
    float numeratorBeta = dxBeta * (xmin + 0.5) + dyBeta * (ymin + 0.5) + x0 * y2 - x2 * y0;
    float numeratorGamma = dxGamma * (xmin + 0.5) + dyGamma * (ymin + 0.5) + x0 * y1 - x1 * y0;
    float denominatorAlpha = dxAlpha * x0 + dyAlpha * y0 + x1 * y2 - x2 * y1;
    float denominatorBeta = dxBeta * x1 + dyBeta * y1 + x0 * y2 - x2 * y0;
    float denominatorGamma = dxGamma * x2 + dyGamma * y2 + x0 * y1 - x1 * y0;
    for (int x = xmin; x <= xmax; ++x) {
        float numeratorAlphaTmp = numeratorAlpha;
        float numeratorBetaTmp = numeratorBeta;
        float numeratorGammaTmp = numeratorGamma;
        for (int y = ymin; y <= ymax; ++y) {
            float alpha = numeratorAlphaTmp / denominatorAlpha;
            float beta = numeratorBetaTmp / denominatorBeta;
            float gamma = numeratorGammaTmp / denominatorGamma;
            if (alpha >= 0 && beta >= 0 && gamma >= 0)
#ifdef SSAA
                render_point(x, y, color, &supersample_target[0], target_w * sample_rate);
#elif defined(MLAA)
                render_point(x, y, color, &mlaa_sample_target[0], target_w);
#else
                render_point(x, y, color, render_target, target_w);
#endif
            numeratorAlphaTmp += dyAlpha;
            numeratorBetaTmp += dyBeta;
            numeratorGammaTmp += dyGamma;
        }
        numeratorAlpha += dxAlpha;
        numeratorBeta += dxBeta;
        numeratorGamma += dxGamma;
    }
}

void SoftwareRendererImp::rasterize_image(float x0, float y0, float x1,
    float y1, Texture& tex)
{
    // Task 6:
    // Implement image rasterization
}

#ifdef SSAA
// resolve samples to render target
Color SoftwareRendererImp::resolve_point(int x, int y)
{
    Color color { 0, 0, 0, 0 };
    int xSample = x * sample_rate;
    int ySample = y * sample_rate;
    int wSample = target_w * sample_rate;
    for (int i = 0; i < sample_rate; i++)
        for (int j = 0; j < sample_rate; j++)
        {
            int index = 4 * (xSample + i + (ySample + j) * wSample);
            color.r += supersample_target[index] / 255.f;
            color.g += supersample_target[index + 1] / 255.f;
            color.b += supersample_target[index + 2] / 255.f;
            color.a += supersample_target[index + 3] / 255.f;
        }
    return color * (1.f / (sample_rate * sample_rate));
}
#endif

#ifdef MLAA
enum class PatternKind : unsigned char {
    NONE,
    H,
    B,
    T,
    L
};
struct PatternInfo {
    int x;
    int y;
    int distance;
    PatternKind start;
    PatternKind end;
};
struct PatternBias {
    float from;
    float to;
};
union AreaInfo {
    struct {
        float top;
        float down;
    };
    struct {
        float left;
        float right;
    };
};
PatternBias analyse_pattern(PatternKind start, PatternKind end)
{
    if (start == PatternKind::H) {
        if (end == PatternKind::H)
            return { 0, 0 };
        else if (end == PatternKind::B)
            return { 0.5, -0.5 };
        else if (end == PatternKind::T)
            return { -0.5, 0.5 };
        else
            return { 0, 0 };
    } else if (start == PatternKind::B) {
        if (end == PatternKind::H)
            return { -0.5, 0.5 };
        else if (end == PatternKind::B)
            return { -0.5, -0.5 };
        else if (end == PatternKind::T)
            return { -0.5, 0.5 };
        else
            return { -0.5, 0 };
    } else if (start == PatternKind::T) {
        if (end == PatternKind::H)
            return { 0.5, -0.5 };
        else if (end == PatternKind::B)
            return { 0.5, -0.5 };
        else if (end == PatternKind::T)
            return { 0.5, 0.5 };
        else
            return { 0.5, 0 };
    } else {
        if (end == PatternKind::H)
            return { 0, 0 };
        else if (end == PatternKind::B)
            return { 0, -0.5 };
        else if (end == PatternKind::T)
            return { 0, 0.5 };
        else
            return { 0, 0 };
    }
}
void caculate_area_of_pattern(int distance, PatternBias& bias, std::vector<AreaInfo>& areaInfo)
{
    float height = 0;
    float bottom = 0;
    if (bias.from == 0 && bias.to == 0) // O pattern -- LL Lh HL HH pattern
        return;
    else if (bias.from == 0) {
        height = bias.to;
        bottom = distance;
    } else if (bias.to == 0) {
        height = bias.from;
        bottom = distance;
    } else {
        height = bias.from;
        bottom = distance / 2.0f;
    }
    float triArea = std::abs(height) * bottom / 2;

    if (bias.from == 0) { // LB LT pattern
        for (int i = 0; i < distance; ++i) {
            float ratio = ((i + 1) / bottom) * ((i + 1) / bottom) - (i / bottom) * (i / bottom);
            float area = bias.to * 2 * triArea * ratio;
            if (area > 0)
                areaInfo.push_back({ area, 0 });
            else
                areaInfo.push_back({ 0, -area });
        }
    } else if (bias.to == 0) { // BL TL pattern
        for (int i = 0; i < distance; ++i) {
            float ratio = ((bottom - i) / bottom) * ((bottom - i) / bottom) - ((bottom - i - 1) / bottom) * ((bottom - i - 1) / bottom);
            float area = bias.from * 2 * (triArea * ratio);
            if (area > 0)
                areaInfo.push_back({ area, 0 });
            else
                areaInfo.push_back({ 0, -area });
        }
    } else if (distance % 2 != 0) { // bottom line length is not integer
        int boundary = std::floor(bottom);
        for (int i = 0; i < distance; ++i) {
            if (i == boundary)
            {
                if (bias.from != bias.to) // Z pattern
                    areaInfo.push_back({ 0,0 });
                else{ // U pattern
                    float ratio = (0.5 / bottom) * (0.5 / bottom);
                    float area = bias.from * 2 * triArea * ratio * 2;
                    if (area > 0)
                        areaInfo.push_back({ area, 0 });
                    else
                        areaInfo.push_back({ 0, -area });
                }
            }
            else if (i < boundary) { // U and Z degrade into BL TL pattern
                float ratio = ((bottom - i) / bottom) * ((bottom - i) / bottom) - ((bottom - i - 1) / bottom) * ((bottom - i - 1) / bottom);
                float area = bias.from * 2 * triArea * ratio;
                if (area > 0)
                    areaInfo.push_back({ area, 0 });
                else
                    areaInfo.push_back({ 0, -area });
            }
            else { // U and Z pattern
                float ratio = ((i - boundary + 0.5) / bottom) * ((i - boundary + 0.5) / bottom) - ((i - boundary - 0.5) / bottom) * ((i - boundary - 0.5) / bottom);
                float area = bias.to * 2 * triArea * ratio;
                if (area > 0)
                    areaInfo.push_back({ area, 0 });
                else
                    areaInfo.push_back({ 0, -area });
            }
        }
        
    } else { // bottom line length is integer
        for (int i = 0; i < distance; ++i) {
            if (i < bottom) { //  U and Z degrade into BL TL pattern
                float ratio = ((bottom - i) / bottom) * ((bottom - i) / bottom) - ((bottom - i - 1) / bottom) * ((bottom - i - 1) / bottom);
                float area = bias.from * 2 * triArea * ratio;
                if (area > 0)
                    areaInfo.push_back({ area, 0 });
                else
                    areaInfo.push_back({ 0, -area });
            }
            else { // U and Z pattern
                float ratio = ((i - bottom + 1) / bottom) * ((i - bottom + 1) / bottom) - ((i - bottom) / bottom) * ((i - bottom) / bottom);
                float area = bias.to * 2 * triArea * ratio;
                if (area > 0)
                    areaInfo.push_back({ area, 0 });
                else
                    areaInfo.push_back({ 0, -area });
            }
        }
    }
}
#endif // MLAA

void SoftwareRendererImp::resolve(void)
{

    // Implement supersampling
    // You may also need to modify other functions marked with "Task 4".
#ifdef SSAA
    for (size_t y = 0; y < target_h; ++y) {
        for (size_t x = 0; x < target_w; ++x) {
            int index = 4 * (x + y * target_w);
            Color color = resolve_point(x, y); // average color within supersample_target sample_rate * sample_rate area
            if (color == Color::White)
                color = { render_target_tmp[index] / 255.f,
                    render_target_tmp[index + 1] / 255.f,
                    render_target_tmp[index + 2] / 255.f,
                    render_target_tmp[index + 3] / 255.f };
            render_point(x, y, color, render_target, target_w);
        }
    }
#elif defined(MLAA)
    // decide if the pixel compared with another pixel has edge
    static const float threshold = 0.1 * 255;   

    // detect if enable mlaa
    if (!mlaa)
    {
        for (int y = 0; y < target_h; ++y) {
            for (int x = 0; x < target_w; ++x) {
                int index = x + y * target_w;
                Color color{ mlaa_sample_target[4 * index] / 255.f, mlaa_sample_target[4 * index + 1] / 255.f, mlaa_sample_target[4 * index + 2] / 255.f, mlaa_sample_target[4 * index + 3] / 255.f };
                if (color == Color::White)
                    color = { render_target_tmp[4 * index] / 255.f,
                    render_target_tmp[4 * index + 1] / 255.f,
                    render_target_tmp[4 * index + 2] / 255.f,
                    render_target_tmp[4 * index + 3] / 255.f };
                render_point(x, y, color, render_target, target_w);
            }
        }
        return;
    }

    // turn color image into gray image by formula Y = 0.2126 * R + 0.7152 * G + 0.0722 * B
    std::vector<unsigned char> gray_sample_target;
    gray_sample_target.resize(target_w * target_h);
    for (int y = 0; y < target_h; ++y) {
        for (int x = 0; x < target_w; ++x) {
            int index = x + y * target_w;
            gray_sample_target[index] = 0.2126f * mlaa_sample_target[4 * index] + 
                                        0.7152f * mlaa_sample_target[4 * index + 1] + 
                                        0.0722f * mlaa_sample_target[4 * index + 2];
        }
    }

    // find edges
    std::vector<unsigned char> edge_buffer;
    edge_buffer.resize(3 * target_w * target_h);
    memset(&edge_buffer[0], 0, 3 * target_w * target_h);
    for (int y = 0; y < target_h; ++y) {
        int yAxi = y * target_w;
        for (int x = 1; x < target_w; ++x) {
            int index = x + yAxi;
            //  compare the current pixel with left one and if distance greater than threshold we set current pixel R channel to 1
            if (std::abs(gray_sample_target[index] - gray_sample_target[index - 1]) > threshold)
                edge_buffer[3 * index] = 255;
        }
    }
    for (int y = 1; y < target_h; ++y) {
        int yAxi = y * target_w;
        for (int x = 0; x < target_w; ++x) {
            int index = x + yAxi;
            // compare the current pixel with top one and if distance greater than threshold we set current pixel G channel to 1
            if (std::abs(gray_sample_target[index] - gray_sample_target[index - target_w]) > threshold)
                edge_buffer[3 * index + 1] = 255;
        }
    }

    std::vector<PatternInfo> patternsX;
    for (int y = 1; y < target_h; ++y) {
        int yAxi = y * target_w;
        for (int x = 0; x < target_w;) {
            int index = x + yAxi;
            //  X-axis offset from the current pixel to the next pixel to query
            int bias = 1;

            if (edge_buffer[3 * index + 1]) {
                PatternInfo pattern;
                pattern.x = x;
                pattern.y = y;

                // detect start pattern
                auto up = edge_buffer[3 * (index - target_w)];
                auto down = edge_buffer[3 * index];
                if (down && up)
                    pattern.start = PatternKind::H;
                else if (up)
                    pattern.start = PatternKind::B;
                else if (down)
                    pattern.start = PatternKind::T;
                else
                    pattern.start = PatternKind::L;

                // detect base line length and end pattern
                pattern.end = PatternKind::NONE;
                for (; bias < target_w - x; ++bias) {
                    // continue if no H or B or T and has G edge
                    auto up = edge_buffer[3 * (index - target_w + bias)];
                    auto down = edge_buffer[3 * (index + bias)];
                    if (down && up) {
                        pattern.end = PatternKind::H;
                        break;
                    }
                    else if (up) {
                        pattern.end = PatternKind::B;
                        break;
                    }
                    else if (down) {
                        pattern.end = PatternKind::T;
                        break;
                    }
                    else if (!edge_buffer[3 * (index + bias) + 1])
                    {
                        pattern.end = PatternKind::L;
                        break;
                    }
                }
                if (pattern.end == PatternKind::NONE)
                    pattern.end = PatternKind::L;
                pattern.distance = bias;

                patternsX.push_back(pattern);
            }
            x += bias;
        }
    }
    std::vector<PatternInfo> patternsY;
    for (int x = 1; x < target_w; ++x) {
        for (int y = 0; y < target_h;) {
            int index = x + y * target_w;
            //  Y-axis offset from the current pixel to the next pixel to query
            int bias = 1;

            if (edge_buffer[3 * index]) {
                PatternInfo pattern;
                pattern.x = x;
                pattern.y = y;

                // detect start pattern
                auto left = edge_buffer[3 * (index - 1) + 1];
                auto right = edge_buffer[3 * index + 1];
                if (left && right)
                    pattern.start = PatternKind::H;
                else if (left)
                    pattern.start = PatternKind::B;
                else if (right)
                    pattern.start = PatternKind::T;
                else
                    pattern.start = PatternKind::L;

                // detect base line length and end pattern
                pattern.end = PatternKind::NONE;
                for (; bias < target_h - y; ++bias) {
                    // continue if no H or B or T and has R edge
                    auto left = edge_buffer[3 * (index - 1 + bias * target_w) + 1];
                    auto right = edge_buffer[3 * (index + bias * target_w) + 1];
                    if (left && right) {
                        pattern.end = PatternKind::H;
                        break;
                    } else if (left) {
                        pattern.end = PatternKind::B;
                        break;
                    } else if (right) {
                        pattern.end = PatternKind::T;
                        break;
                    } else if (!edge_buffer[3 * (index + bias * target_w)])
                    {
                        pattern.end = PatternKind::L;
                        break;
                    }
                }
                if (pattern.end == PatternKind::NONE)
                    pattern.end = PatternKind::L;
                pattern.distance = bias;

                patternsY.push_back(pattern);
            }
            y += bias;
        }
    }
    
    std::vector<unsigned char> weight_buffer;
    weight_buffer.resize(4 * target_w * target_h);
    memset(&weight_buffer[0], 0, 4 * target_w * target_h);
    for (auto& patternInfo : patternsX) {
        std::vector<AreaInfo> areaInfos;
        caculate_area_of_pattern(patternInfo.distance, analyse_pattern(patternInfo.start, patternInfo.end), areaInfos);
        int size = areaInfos.size();
        for (int i = 0; i < size; ++i) {
            weight_buffer[4 * (patternInfo.x + i + patternInfo.y * target_w)] = areaInfos[i].top * 255;
            weight_buffer[4 * (patternInfo.x + i + patternInfo.y * target_w) + 1] = areaInfos[i].down * 255;
        }
    }
    for (auto& patternInfo : patternsY) {
        std::vector<AreaInfo> areaInfos;
        caculate_area_of_pattern(patternInfo.distance, analyse_pattern(patternInfo.start, patternInfo.end), areaInfos);
        int size = areaInfos.size();
        for (int i = 0; i < size; ++i) {
            weight_buffer[4 * (patternInfo.x + (patternInfo.y + i) * target_w) + 2] = areaInfos[i].left * 255;
            weight_buffer[4 * (patternInfo.x + (patternInfo.y + i) * target_w) + 3] = areaInfos[i].right * 255;
        }
    }

    // loop for color blend and show final mlaa image
    // current pixel color is blended by formual ( Ctop * wTop + Cdown * wDown + Cleft * wLeft + Cright * wRight ) / 2, wRop and wDown 
    for (size_t y = 1; y < target_h-1; ++y) {
        int yAxi = y * target_w;
        for (size_t x = 1; x < target_w-1; ++x) {
            int index = x + yAxi;
            Color old = { mlaa_sample_target[4 * index] / 255.f, mlaa_sample_target[4 * index + 1] / 255.f, mlaa_sample_target[4 * index + 2] / 255.f, mlaa_sample_target[4 * index + 3] / 255.f };
            Color top = { mlaa_sample_target[4 * (index - target_w)] / 255.f, mlaa_sample_target[4 * (index - target_w) + 1] / 255.f, mlaa_sample_target[4 * (index - target_w) + 2] / 255.f, mlaa_sample_target[4 * (index - target_w) + 3] / 255.f };
            Color down = { mlaa_sample_target[4 * (index + target_w)] / 255.f, mlaa_sample_target[4 * (index + target_w) + 1] / 255.f, mlaa_sample_target[4 * (index + target_w) + 2] / 255.f, mlaa_sample_target[4 * (index + target_w) + 3] / 255.f };
            Color left = { mlaa_sample_target[4 * (index - 1)] / 255.f, mlaa_sample_target[4 * (index - 1) + 1] / 255.f, mlaa_sample_target[4 * (index - 1) + 2] / 255.f, mlaa_sample_target[4 * (index - 1) + 3] / 255.f };
            Color right = { mlaa_sample_target[4 * (index + 1)] / 255.f, mlaa_sample_target[4 * (index + 1) + 1] / 255.f, mlaa_sample_target[4 * (index + 1) + 2] / 255.f,mlaa_sample_target[4 * (index + 1) + 3] / 255.f };
            float wTop = weight_buffer[4 * index] / 255.f;
            float wLeft = weight_buffer[4 * index + 2] / 255.f;
            float wDown = weight_buffer[4 * (index + target_w) + 1] / 255.f;
            float wRight = weight_buffer[4 * (index + 1) + 3] / 255.f;
            Color newX = (1 - wTop - wDown) * old + wTop * top + wDown * down;
            Color newY = (1 - wLeft - wRight) * old + wLeft * left + wRight * right;
            Color newColor = (newX + newY) * 0.5;
            if (newColor == Color::White)
                newColor = { render_target_tmp[4 * index] / 255.f,
                    render_target_tmp[4 * index + 1] / 255.f,
                    render_target_tmp[4 * index + 2] / 255.f,
                    render_target_tmp[4 * index + 3] / 255.f };
            render_point(x, y, newColor, render_target, target_w);
        }
    }

    // loop for showing gray image
    /*for (size_t y = 0; y < target_h; y++)
    {
        int yAxi = y * target_w;
        for (size_t x = 0; x < target_w; ++x) {
            int index = x + yAxi;
            Color color{ gray_sample_target[index] / 255.f, gray_sample_target[index] / 255.f, gray_sample_target[index] / 255.f };
            if (color == Color::White)
                color = { render_target_tmp[4 * index] / 255.f,
                          render_target_tmp[4 * index + 1] / 255.f,
                          render_target_tmp[4 * index + 2] / 255.f,
                          render_target_tmp[4 * index + 3] / 255.f };
            render_point(x, y, color, render_target, target_w);
        }
    }*/

     // loop for showing edge image
     /*for (size_t y = 0; y < target_h; y++) {
         int yAxi = y * target_w;
         for (size_t x = 0; x < target_w; ++x) {
             int index = x + yAxi;
             Color color { edge_buffer[3 * index] / 255.f, edge_buffer[3 * index + 1] / 255.f, edge_buffer[3 * index + 2] / 255.f };
             if (color == Color::Black)
                 color = { render_target_tmp[4 * index] / 255.f,
                     render_target_tmp[4 * index + 1] / 255.f,
                     render_target_tmp[4 * index + 2] / 255.f,
                     render_target_tmp[4 * index + 3] / 255.f };
             render_point(x, y, color, render_target, target_w);
         }
     }*/

     // loop for showing weight image
     /*Color defaultColor{ 0,0,0,1 };
     for (size_t y = 0; y < target_h; y++) {
         int yAxi = y * target_w;
         for (size_t x = 0; x < target_w; ++x) {
             int index = x + yAxi;
             Color color{ weight_buffer[4 * index] / 255.f, weight_buffer[4 * index + 1] / 255.f, weight_buffer[4 * index + 2] / 255.f, weight_buffer[4 * index + 3] / 255.f };
             color.a = 1;
             if (color == defaultColor)
                 color = { render_target_tmp[4 * index] / 255.f,
                     render_target_tmp[4 * index + 1] / 255.f,
                     render_target_tmp[4 * index + 2] / 255.f,
                     render_target_tmp[4 * index + 3] / 255.f };
             render_point(x, y, color, render_target, target_w);
         }
     }*/

#endif
    return;
}

} // namespace CMU462
