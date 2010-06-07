/*
 * Copyright (c) 2008 The Jackson Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.jax.haplotype.analysis.visualization;

import org.jax.geneticutil.data.IndexedSnpInterval;

/**
 * A class for defining a visualization interval that allows for some of the
 * indices to be highlighted
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HighlightedSnpInterval extends IndexedSnpInterval
{
    private final int[] indicesToHighlight;

    /**
     * Constructor
     * @param startIndex
     *          see {@link #getExtentInIndices()}
     * @param extentInIndices
     *          see {@link #getExtentInIndices()}
     * @param indicesToHighlight
     *          the indices to highlight
     */
    public HighlightedSnpInterval(
            int startIndex,
            int extentInIndices,
            int[] indicesToHighlight)
    {
        super(startIndex, extentInIndices);
        this.indicesToHighlight = indicesToHighlight;
    }
    
    /**
     * Getter for the index to highlight
     * @return the indicesToHighlight
     */
    public int[] getIndicesToHighlight()
    {
        return this.indicesToHighlight;
    }
}
