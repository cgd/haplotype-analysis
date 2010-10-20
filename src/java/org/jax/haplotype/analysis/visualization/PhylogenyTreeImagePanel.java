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

import java.awt.Graphics;
import java.awt.image.BufferedImage;

import javax.swing.JPanel;

import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;

/**
 * A panel for rendering phylogeny trees
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class PhylogenyTreeImagePanel extends JPanel
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 6701039578229770767L;
    
    private final PhylogenyTreeImageFactory treeImageFactory;
    
    private volatile PhylogenyTreeNode phylogenyTree;
    
    

    /**
     * Constructor
     * @param treeImageFactory
     *          the tree image factory
     */
    public PhylogenyTreeImagePanel(PhylogenyTreeImageFactory treeImageFactory)
    {
        this(treeImageFactory, null);
    }



    /**
     * Constructor
     * @param treeImageFactory
     *          the phylogeny tree image factory to use for rendering
     * @param phylogenyTree
     *          the phylogeny tree node to render
     */
    public PhylogenyTreeImagePanel(
            PhylogenyTreeImageFactory treeImageFactory,
            PhylogenyTreeNode phylogenyTree)
    {
        this.treeImageFactory = treeImageFactory;
        this.phylogenyTree = phylogenyTree;
    }
    
    /**
     * Setter for the phylogeny tree to render
     * @param phylogenyTree the phylogenyTree to set
     */
    public void setPhylogenyTree(PhylogenyTreeNode phylogenyTree)
    {
        this.phylogenyTree = phylogenyTree;
        this.repaint();
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void paintComponent(Graphics g)
    {
        super.paintComponent(g);
        PhylogenyTreeNode phylogenyTree = this.phylogenyTree;
        if(phylogenyTree != null)
        {
            BufferedImage treeImage = this.treeImageFactory.createPhylogenyImage(
                    phylogenyTree,
                    this.getWidth(),
                    this.getHeight());
            g.drawImage(treeImage, 0, 0, this);
        }
    }
}
