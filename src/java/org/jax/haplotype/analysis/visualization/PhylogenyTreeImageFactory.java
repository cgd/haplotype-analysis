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

import java.awt.image.BufferedImage;

import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;

/**
 * Interface for a class that can generate a phylogeny image.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface PhylogenyTreeImageFactory
{
    /**
     * Create a phylogeny image
     * @param phylogenyTree
     *          the tree to create an image from
     * @param imageWidth
     *          the image width
     * @param imageHeight
     *          the image height
     * @return
     *          the image
     */
    public BufferedImage createPhylogenyImage(
            PhylogenyTreeNode phylogenyTree,
            int imageWidth,
            int imageHeight);
}
