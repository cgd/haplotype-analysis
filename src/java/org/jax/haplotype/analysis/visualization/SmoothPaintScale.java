/*
 * Copyright (c) 2010 The Jackson Laboratory
 * 
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.jax.haplotype.analysis.visualization;

import java.awt.Color;
import java.io.Serializable;

import org.jfree.chart.renderer.PaintScale;

/**
 * A paint scale that can be used for smooth color transitions
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class SmoothPaintScale implements PaintScale, Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -9128634395046425457L;

    private final Color upperSmoothColor;
    
    private final Color lowerSmoothColor;
    
    private final double upperBound;
    
    private final double lowerBound;
    
    /**
     * Constructor
     * @param upperBound
     *          the upper bound value
     * @param lowerBound
     *          the lower bound value
     * @param lowerSmoothColor
     *          the color to smooth toward as we reach the lower part of the
     *          range
     * @param upperSmoothColor
     *          the color to smooth toward as we reach the upper part of the
     *          range
     */
    public SmoothPaintScale(
            double lowerBound,
            double upperBound,
            Color lowerSmoothColor,
            Color upperSmoothColor)
    {
        this.upperBound = upperBound;
        this.lowerBound = lowerBound;
        this.lowerSmoothColor = lowerSmoothColor;
        this.upperSmoothColor = upperSmoothColor;
    }

    /**
     * {@inheritDoc}
     */
    public double getLowerBound()
    {
        return this.lowerBound;
    }

    /**
     * {@inheritDoc}
     */
    public Color getPaint(double value)
    {
        double totalDiff = this.upperBound - this.lowerBound;
        
        double upperColorWeight = (value - this.lowerBound) / totalDiff;
        double lowerColorWeight = 1.0 - upperColorWeight;
        
        int red = this.boundColorComponent(
            this.upperSmoothColor.getRed() * upperColorWeight +
            this.lowerSmoothColor.getRed() * lowerColorWeight);
        int green = this.boundColorComponent(
            this.upperSmoothColor.getGreen() * upperColorWeight +
            this.lowerSmoothColor.getGreen() * lowerColorWeight);
        int blue = this.boundColorComponent(
            this.upperSmoothColor.getBlue() * upperColorWeight +
            this.lowerSmoothColor.getBlue() * lowerColorWeight);
        
        return new Color(red, green, blue);
    }
    
    /**
     * make sure that we have a valid color even if the value used to
     * create the color goes beyond the specified lower and upper bounds
     * @param colorComponent
     *          the real-valued color component (on the 0-255 scale)
     * @return
     *          an int version of colorComponent which is bound between
     *          0 and 255
     */
    private int boundColorComponent(double colorComponent)
    {
        int roundedValue = (int)Math.round(colorComponent);
        
        if(roundedValue < 0)
        {
            return 0;
        }
        else if(roundedValue > 255)
        {
            return 255;
        }
        else
        {
            return roundedValue;
        }
    }

    /**
     * {@inheritDoc}
     */
    public double getUpperBound()
    {
        return this.upperBound;
    }
}
