#ifndef _PLOTUTILS_H_
#define _PLOTUTILS_H_

#include <vector>
#include <string>

namespace PlotUtils
{
    void SetStyle();
    void myText( double x, double y, int color, const char * text, const float size = 0.03 );

} // namespace PlotUtils

#endif // _PLOTUTILS_H_