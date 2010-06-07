/*
 * Copyright (c) 2008 The Jackson Laboratory
 *
 * Permission is hereby granted, free of charge, to any person obtaining  a copy
 * of this software and associated documentation files (the  "Software"), to
 * deal in the Software without restriction, including  without limitation the
 * rights to use, copy, modify, merge, publish,  distribute, sublicense, and/or
 * sell copies of the Software, and to  permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be  included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF  MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE  SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
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
