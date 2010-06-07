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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Paint;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.font.FontRenderContext;
import java.awt.geom.AffineTransform;
import java.awt.geom.Area;
import java.awt.geom.Dimension2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.RoundRectangle2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.jax.haplotype.phylogeny.data.PhylogenyTreeEdge;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeEdgeWithRealValue;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;
import org.jax.util.TextWrapper;
import org.jax.util.gui.Dimension2DDouble;
import org.jax.util.gui.Java2DUtils;
import org.jfree.chart.renderer.PaintScale;

/**
 * A phylogeny tree image factory that doesn't depend on external
 * tree rendering libraries
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class SimplePhylogenyTreeImageFactory implements PhylogenyTreeImageFactory
{
    private static final Logger LOG = Logger.getLogger(
            SimplePhylogenyTreeImageFactory.class.getName());
    
    private static final Font LABEL_FONT = new Font("SansSerif", Font.BOLD, 12);
    
    private static final int LABEL_BORDER_ROUNDING_ARC_SIZE = 10;
    
    private static final double LEAF_ANGLE_WEIGHTING = 2;
    
    // TODO imlement some kind of strain weighting scheme
    private static final double STRAIN_ANGLE_WEIGHTING = 0;
    
    private static final Color FOREGROUND_COLOR = Color.BLACK;
    
    private static final Color SHADOW_COLOR = Color.LIGHT_GRAY;
    
    private static final Color BACKGROUND_COLOR =
        new Color(1.0F, 1.0F, 1.0F, 0.5F);
    
    private static final Stroke LINE_STROKE =
        new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
    
    private static final AffineTransform SHADOW_TRANSFORM =
        new AffineTransform();
    static
    {
        // this is the transform that we use to draw the pretty shadow
        // background
        SHADOW_TRANSFORM.translate(3.0, 3.0);
    }
    
    private static final int BORDER_WHITE_SPACE = 10;

    private static final int NODE_LABEL_WRAP_LIMIT = 25;
    
    private final PaintScale paintScale;
    
    /**
     * Constructor
     */
    public SimplePhylogenyTreeImageFactory()
    {
        this(null);
    }
    
    /**
     * Constructor
     * @param paintScale
     *          the paint scale to use
     */
    public SimplePhylogenyTreeImageFactory(PaintScale paintScale)
    {
        this.paintScale = paintScale;
    }

    /**
     * {@inheritDoc}
     */
    public BufferedImage createPhylogenyImage(
            PhylogenyTreeNode phylogenyTree,
            int imageWidth,
            int imageHeight)
    {
        if(LOG.isLoggable(Level.FINE))
        {
            LOG.fine(
                    "Creating phylogeny image for: " +
                    phylogenyTree.resolveToSingleStrainLeafNodes(0.0).toNewickFormat());
        }
        
        BufferedImage bi = new BufferedImage(
                imageWidth,
                imageHeight,
                BufferedImage.TYPE_4BYTE_ABGR);
        Graphics2D graphics = bi.createGraphics();
        graphics.setStroke(LINE_STROKE);
        graphics.setRenderingHint(
                RenderingHints.KEY_RENDERING,
                RenderingHints.VALUE_RENDER_QUALITY);
        graphics.setRenderingHint(
                RenderingHints.KEY_ANTIALIASING,
                RenderingHints.VALUE_ANTIALIAS_ON);
        graphics.setColor(Color.BLACK);
        
        this.paintPhylogenyTree(
                graphics,
                phylogenyTree,
                imageWidth,
                imageHeight);
        
        return bi;
    }
    
    /**
     * Paint the given tree
     * @param graphics
     *          the graphics to paint with
     * @param phylogenyTree
     *          the tree to paint
     */
    private void paintPhylogenyTree(
            Graphics2D graphics,
            PhylogenyTreeNode phylogenyTree,
            int imageWidth,
            int imageHeight)
    {
        VisualTreeNode treeLayout = this.createTreeLayout(phylogenyTree);
        this.transformTreeLayout(
                treeLayout,
                imageWidth,
                imageHeight,
                graphics.getFontRenderContext());
        this.paintPhylogenyTree(
                graphics,
                treeLayout);
    }

    /**
     * Paint the given tree layout
     * @param graphics
     *          the graphics to paint with
     * @param treeLayout
     *          the layout to paint
     */
    private void paintPhylogenyTree(
            Graphics2D graphics,
            VisualTreeNode treeLayout)
    {
        int childNodeCount = treeLayout.getChildNodes().size();
        for(int i = 0; i < childNodeCount; i++)
        {
            VisualTreeNode visualChild = treeLayout.getChildNodes().get(i);
            
            Shape branchShape = new Line2D.Double(
                    treeLayout.getPosition(),
                    visualChild.getPosition());
            Shape branchShadowShape = SHADOW_TRANSFORM.createTransformedShape(
                    branchShape);
            
            graphics.setColor(SHADOW_COLOR);
            graphics.draw(branchShadowShape);
            
            if(this.paintScale != null)
            {
                PhylogenyTreeEdge phylogenyEdge =
                    treeLayout.getPhylogenyTreeNode().getChildEdges().get(i);
                if(phylogenyEdge instanceof PhylogenyTreeEdgeWithRealValue)
                {
                    PhylogenyTreeEdgeWithRealValue phylogenyEdgeWithValue =
                        (PhylogenyTreeEdgeWithRealValue)phylogenyEdge;
                    Paint paint = this.paintScale.getPaint(
                            phylogenyEdgeWithValue.getRealValue());
                    graphics.setPaint(paint);
                }
                else
                {
                    graphics.setColor(FOREGROUND_COLOR);
                }
            }
            else
            {
                graphics.setColor(FOREGROUND_COLOR);
            }
            graphics.draw(branchShape);
            
            // recurse
            this.paintPhylogenyTree(
                    graphics,
                    visualChild);
        }
        
        if(!treeLayout.getPhylogenyTreeNode().getStrains().isEmpty())
        {
            Shape textShape =
                this.getLabelShape(treeLayout, graphics.getFontRenderContext());
            Shape borderShape = this.getLabelBorder(textShape);
            
            graphics.setColor(BACKGROUND_COLOR);
            graphics.fill(borderShape);
            
            graphics.setColor(SHADOW_COLOR);
            Area borderShadowShape =
                new Area(SHADOW_TRANSFORM.createTransformedShape(borderShape));
            borderShadowShape.subtract(new Area(borderShape));
            graphics.draw(borderShadowShape);
            
            graphics.setColor(FOREGROUND_COLOR);
            graphics.draw(borderShape);
            
            graphics.fill(textShape);
        }
    }
    
    /**
     * Get the label shape for the given tree node
     * @param treeNode
     *          the tree node
     * @param frc
     *          the font rendering context
     * @return
     *          the label shape
     */
    private Shape getLabelShape(VisualTreeNode treeNode, FontRenderContext frc)
    {
        // convert the strain list to a comma separated list then wrap it
        List<String> strainList = treeNode.getPhylogenyTreeNode().getStrains();
        int strainCount = strainList.size();
        StringBuffer commaSeparatedStrains = new StringBuffer();
        for(int i = 0; i < strainCount; i++)
        {
            if(i >= 1)
            {
                commaSeparatedStrains.append(", ");
            }
            commaSeparatedStrains.append(strainList.get(i));
        }
        String[] wrappedStrains = TextWrapper.wrapText(
                commaSeparatedStrains.toString(),
                NODE_LABEL_WRAP_LIMIT);
        
        Shape textShape = Java2DUtils.createCenteredMultilineTextShape(
                wrappedStrains,
                LABEL_FONT,
                frc);
        AffineTransform transform = new AffineTransform();
        transform.translate(
                treeNode.getPosition().getX(),
                treeNode.getPosition().getY());
        textShape = transform.createTransformedShape(textShape);
        
        return textShape;
    }
    
    /**
     * Get the border shape for the given label shape
     * @param labelShape
     *          the label shape that we're going to draw a border around
     * @return
     *          the border shape
     */
    private Shape getLabelBorder(Shape labelShape)
    {
        Rectangle2D labelBounds = labelShape.getBounds2D();
        
        return new RoundRectangle2D.Double(
                labelBounds.getX() - LABEL_BORDER_ROUNDING_ARC_SIZE,
                labelBounds.getY() - LABEL_BORDER_ROUNDING_ARC_SIZE,
                labelBounds.getWidth() + (LABEL_BORDER_ROUNDING_ARC_SIZE * 2),
                labelBounds.getHeight() + (LABEL_BORDER_ROUNDING_ARC_SIZE * 2),
                LABEL_BORDER_ROUNDING_ARC_SIZE,
                LABEL_BORDER_ROUNDING_ARC_SIZE);
    }
    
    /**
     * Transform the tree layout so that it fits nicely in the given image
     * dimensions
     * @param treeLayout
     *          the layout to transform
     * @param imageWidth
     *          the image width
     * @param imageHeight
     *          the image height
     */
    private void transformTreeLayout(
            VisualTreeNode treeLayout,
            int imageWidth,
            int imageHeight,
            FontRenderContext frc)
    {
        Dimension2D maximalNodeLabelDimension = this.calculateMaximalNodeDimension(
                treeLayout,
                frc);
        double widthBuffer = maximalNodeLabelDimension.getWidth() + BORDER_WHITE_SPACE;
        double heightBuffer = maximalNodeLabelDimension.getHeight() + BORDER_WHITE_SPACE;
        
        // perform rotation to improve the use of space
        {
            // center around 0, 0
            VisualTreeNode[] mostDistantPair =
                this.getMostDistantNodePair(treeLayout);
            Point2D distantPoint1 = mostDistantPair[0].getPosition();
            Point2D distantPoint2 = mostDistantPair[1].getPosition();
            double xDiff = distantPoint1.getX() - distantPoint2.getX();
            double yDiff = distantPoint1.getY() - distantPoint2.getY();
            this.translateTreeLayout(
                    treeLayout,
                    (xDiff / 2.0) - distantPoint1.getX(),
                    (yDiff / 2.0) - distantPoint2.getY());
            
            // rotate
            double thetaRadians = Math.atan2(yDiff, xDiff);
            
            if(imageWidth >= imageHeight)
            {
                this.rotateTreeLayout(treeLayout, -thetaRadians);
            }
            else
            {
                this.rotateTreeLayout(treeLayout, (Math.PI / 2.0 - thetaRadians));
            }
        }
        
        Rectangle2D boundingRectangle = this.calculateBounds(treeLayout, null);
        
        // center around the middle of the display area
        this.translateTreeLayout(
                treeLayout,
                -boundingRectangle.getX(),
                -boundingRectangle.getY());
        
        // grow the image to fill a larger area
        double xScale = (imageWidth - widthBuffer) / boundingRectangle.getWidth();
        double yScale = (imageHeight - heightBuffer) / boundingRectangle.getHeight();
        double smallerScale = Math.min(xScale, yScale);
        
        this.scaleTreeLayout(treeLayout, smallerScale);
        
        // center around the middle of the display area
        boundingRectangle = this.calculateBounds(treeLayout, null);
        this.translateTreeLayout(
                treeLayout,
                ((imageWidth - boundingRectangle.getWidth()) / 2.0) - boundingRectangle.getX(),
                ((imageHeight - boundingRectangle.getHeight()) / 2.0) - boundingRectangle.getY());
    }

    /**
     * Find the biggest width and height of tree node labels
     * @param treeLayout
     *          the tree layout
     * @param frc
     *          the {@link FontRenderContext}
     * @return
     *          the biggest width and height
     */
    private Dimension2D calculateMaximalNodeDimension(
            VisualTreeNode treeLayout,
            FontRenderContext frc)
    {
        Dimension2DDouble maximalNodeDimension = new Dimension2DDouble();
        
        this.calculateMaximalNodeDimensionRecursive(
                maximalNodeDimension,
                treeLayout,
                frc);
        
        return maximalNodeDimension;
    }

    /**
     * A recursive version of
     * {@link #calculateMaximalNodeDimension(VisualTreeNode, FontRenderContext)}
     * @param maximalNodeDimension
     *          the biggest width and height
     * @param treeLayout
     *          the tree layout
     * @param frc
     *          the {@link FontRenderContext}
     */
    private void calculateMaximalNodeDimensionRecursive(
            Dimension2DDouble maximalNodeDimension,
            VisualTreeNode treeLayout,
            FontRenderContext frc)
    {
        Shape textShape = this.getLabelShape(treeLayout, frc);
        Shape borderShape = this.getLabelBorder(textShape);
        Rectangle2D borderBounds = borderShape.getBounds2D();
        
        if(borderBounds.getWidth() > maximalNodeDimension.getWidth())
        {
            maximalNodeDimension.setWidth(borderBounds.getWidth());
        }
        
        if(borderBounds.getHeight() > maximalNodeDimension.getHeight())
        {
            maximalNodeDimension.setHeight(borderBounds.getHeight());
        }
        
        for(VisualTreeNode childNode: treeLayout.getChildNodes())
        {
            this.calculateMaximalNodeDimensionRecursive(
                    maximalNodeDimension,
                    childNode,
                    frc);
        }
    }

    /**
     * Rotate the given tree layout
     * @param treeLayout
     *          the layout to rotate
     * @param thetaRadians
     *          the rotation in radians
     */
    private void rotateTreeLayout(
            VisualTreeNode treeLayout,
            double thetaRadians)
    {
        double cosTheta = Math.cos(thetaRadians);
        double sinTheta = Math.sin(thetaRadians);
        this.rotateTreeLayoutRecursive(
                treeLayout,
                cosTheta,
                sinTheta);
    }

    /**
     * Recursively rotate the given tree node
     * @param treeNode
     *          the node to recursively rotate
     * @param cosTheta
     *          the cosine of the rotation angle
     * @param sinTheta
     *          the sin of the rotation angle
     */
    private void rotateTreeLayoutRecursive(
            VisualTreeNode treeNode,
            double cosTheta,
            double sinTheta)
    {
        // rotate this node
        double origX = treeNode.getPosition().getX();
        double origY = treeNode.getPosition().getY();
        double rotX = (origX * cosTheta) - (origY * sinTheta);
        double rotY = (origX * sinTheta) + (origY * cosTheta);
        treeNode.getPosition().setLocation(rotX, rotY);
        
        for(VisualTreeNode childNode: treeNode.getChildNodes())
        {
            // recurse
            this.rotateTreeLayoutRecursive(childNode, cosTheta, sinTheta);
        }
    }

    /**
     * Apply a scaling factor to the tree layout
     * @param treeLayout
     *          the layout to scale
     * @param scale
     *          the scaling factor to apply
     */
    private void scaleTreeLayout(VisualTreeNode treeLayout, double scale)
    {
        treeLayout.getPosition().x *= scale;
        treeLayout.getPosition().y *= scale;
        
        for(VisualTreeNode child: treeLayout.getChildNodes())
        {
            this.scaleTreeLayout(child, scale);
        }
    }

    /**
     * Apply a translation to the tree layout
     * @param treeLayout
     *          the layout to translate
     * @param xTranslation
     *          the x translation
     * @param yTranslation
     *          the y translation
     */
    private void translateTreeLayout(
            VisualTreeNode treeLayout,
            double xTranslation,
            double yTranslation)
    {
        treeLayout.getPosition().x += xTranslation;
        treeLayout.getPosition().y += yTranslation;
        
        for(VisualTreeNode child: treeLayout.getChildNodes())
        {
            this.translateTreeLayout(
                    child,
                    xTranslation,
                    yTranslation);
        }
    }

    /**
     * Method for calculating the minimum bounding rectangle for the given
     * tree layout
     * @param treeLayout
     *          the layout to calculate MBR for
     * @param rectangle
     *          the rectangle up to now (should be null initially)
     * @return
     *          the MBR
     */
    private Rectangle2D.Double calculateBounds(VisualTreeNode treeLayout, Rectangle2D.Double rectangle)
    {
        Point2D position = treeLayout.getPosition();
        if(rectangle == null)
        {
            rectangle = new Rectangle2D.Double(
                    position.getX(), position.getY(), 0.0, 0.0);
        }
        else
        {
            if(position.getX() < rectangle.getMinX())
            {
                double xDiff = rectangle.getMinX() - position.getX();
                rectangle.x -= xDiff;
                rectangle.width += xDiff;
            }
            else if(position.getX() > rectangle.getMaxX())
            {
                double xDiff = position.getX() - rectangle.getMaxX();
                rectangle.width += xDiff;
            }
            
            if(position.getY() < rectangle.getMinY())
            {
                double yDiff = rectangle.getMinY() - position.getY();
                rectangle.y -= yDiff;
                rectangle.height += yDiff;
            }
            else if(position.getY() > rectangle.getMaxY())
            {
                double yDiff = position.getY() - rectangle.getMaxY();
                rectangle.height += yDiff;
            }
        }
        
        for(VisualTreeNode childNode: treeLayout.getChildNodes())
        {
            rectangle = this.calculateBounds(childNode, rectangle);
        }
        
        return rectangle;
    }

    /**
     * Create a layout for the given phylogeny
     * @param phylogenyTree
     *          the tree to create a layout for
     * @return
     *          the layout
     */
    private VisualTreeNode createTreeLayout(
            PhylogenyTreeNode phylogenyTree)
    {
        // how many leaves are there?
        List<PhylogenyTreeNode> leafNodes = phylogenyTree.getAllLeafNodes();
        int leafCount = leafNodes.size();
        int leavesUsed = 0;
        if(phylogenyTree.getChildEdges().size() == 1)
        {
            // this is a leaf if you ignore directionality... so remove
            // one slice of the pie to indicate that we've rendered a leaf
            leafCount++;
            leavesUsed++;
        }
        
        // TODO: implement real strain weighting logic
        // how many strains are there?
        List<String> strains = phylogenyTree.getAllStrains();
        int strainCount = strains.size();
        int strainsUsed = phylogenyTree.getStrains().size();
        
        // figure out the angle weightings between strains and leaves
        double weightedSliceCount =
            (leafCount * LEAF_ANGLE_WEIGHTING) +
            (strainCount * STRAIN_ANGLE_WEIGHTING);
        double weightedRadiansPerSlice =
            (2 * Math.PI) / weightedSliceCount;
        double radiansPerLeaf =
            weightedRadiansPerSlice * LEAF_ANGLE_WEIGHTING;
        double radiansPerStrain =
            weightedRadiansPerSlice * STRAIN_ANGLE_WEIGHTING;
        
        VisualTreeNode root = new VisualTreeNode(phylogenyTree);
        this.createTreeLayoutRecursive(
                root,
                radiansPerLeaf,
                leavesUsed);
        
        return root;
    }

    /**
     * A recursive function for building the tree layout
     * @param parentNode
     *          the parent node that we're filling out
     * @param radiansPerLeaf
     *          the number of radians to use per leaf
     * @param leavesAlreadyUsed
     *          the number of leaves that have already been used by
     *          previous calls
     */
    private void createTreeLayoutRecursive(
            VisualTreeNode parentNode,
            double radiansPerLeaf,
            int leavesAlreadyUsed)
    {
        double parentX = parentNode.getPosition().getX();
        double parentY = parentNode.getPosition().getY();
        
        int cumulativeLeavesUsed = leavesAlreadyUsed;
        for(PhylogenyTreeEdge childEdge: parentNode.getPhylogenyTreeNode().getChildEdges())
        {
            // how many leaves are taken up by this child?
            PhylogenyTreeNode childNode = childEdge.getNode();
            List<PhylogenyTreeNode> childLeafNodes =
                childNode.getAllLeafNodes();
            int currChildLeaves = childLeafNodes.size();
            
            // find the angle to use in placing the child
            double currChildAngleRadians =
                ((currChildLeaves / 2.0) + cumulativeLeavesUsed) * radiansPerLeaf;
            
            // find the x & y positions from the angle and parent position
            // Here's the trig we'll use
            // sin(angle) = y/branch_length
            // y = sin(angle)*branch_length
            // cos(angle) = x/branch_length
            // x = cos(angle)*branch_length
            double branchLength = childEdge.getEdgeLength();
            double yOffset =
                Math.sin(currChildAngleRadians) * branchLength;
            double xOffset =
                Math.cos(currChildAngleRadians) * branchLength;
            
            // create the new child and add it
            VisualTreeNode childLayout = new VisualTreeNode(childNode);
            childLayout.getPosition().setLocation(
                    parentX + xOffset,
                    parentY + yOffset);
            parentNode.getChildNodes().add(childLayout);
            
            // recurse
            this.createTreeLayoutRecursive(
                    childLayout,
                    radiansPerLeaf,
                    cumulativeLeavesUsed);
            cumulativeLeavesUsed += currChildLeaves;
        }
    }
    
    /**
     * Find the tree nodes that are farthest apart
     * @param tree
     *          the tree to search through
     * @return
     *          the most distant pair or null if there are  less than
     *          two nodes
     */
    private VisualTreeNode[] getMostDistantNodePair(VisualTreeNode tree)
    {
        List<VisualTreeNode> leaves = this.getTreeLayoutLeaves(tree);
        
        // TODO there may be a faster way to do this
        VisualTreeNode[] mostDistantPair = null;
        double greatestDistance = -1.0;
        int leafCount = leaves.size();
        for(int i = 0; i < leafCount; i++)
        {
            for(int j = i + 1; j < leafCount; j++)
            {
                double currDistance = leaves.get(i).getPosition().distance(
                        leaves.get(j).getPosition());
                
                if(currDistance > greatestDistance)
                {
                    greatestDistance = currDistance;
                    mostDistantPair = new VisualTreeNode[] {
                            leaves.get(i),
                            leaves.get(j)};
                }
            }
        }
        
        return mostDistantPair;
    }
    
    /**
     * Get all of the leaves in the tree layout (including the tree itself if
     * it is a leaf ignoring directionality
     * @param treeNode
     *          the tree
     * @return
     *          the leaves
     */
    private List<VisualTreeNode> getTreeLayoutLeaves(VisualTreeNode treeNode)
    {
        List<VisualTreeNode> leaves = new ArrayList<VisualTreeNode>();
        
        this.getTreeLayoutLeavesRecursive(treeNode, leaves);
        
        // since directionality doesn't matter we need to check if the
        // root is a leaf too
        if(treeNode.getChildNodes().size() == 1)
        {
            leaves.add(treeNode);
        }
        
        return leaves;
    }

    /**
     * Recursively get leaf nodes
     * @param treeNode
     *          the node to recurse into
     * @param leaves
     *          the leaves to add to
     */
    private void getTreeLayoutLeavesRecursive(
            VisualTreeNode treeNode,
            List<VisualTreeNode> leaves)
    {
        for(VisualTreeNode childNode: treeNode.getChildNodes())
        {
            if(childNode.getChildNodes().isEmpty())
            {
                leaves.add(childNode);
            }
            else
            {
                this.getTreeLayoutLeavesRecursive(
                        childNode,
                        leaves);
            }
        }
    }

    /**
     * An class for positioning tree nodes
     */
    private static final class VisualTreeNode
    {
        private final PhylogenyTreeNode phylogenyTreeNode;
        
        private Point2D.Double position = new Point2D.Double();
        
        private List<VisualTreeNode> childNodes = new ArrayList<VisualTreeNode>();
        
        public VisualTreeNode(PhylogenyTreeNode phylogenyTreeNode)
        {
            this.phylogenyTreeNode = phylogenyTreeNode;
        }

        public List<VisualTreeNode> getChildNodes()
        {
            return this.childNodes;
        }
        
        public Point2D.Double getPosition()
        {
            return this.position;
        }
        
        public void setPosition(Point2D.Double position)
        {
            this.position = position;
        }
        
        public PhylogenyTreeNode getPhylogenyTreeNode()
        {
            return this.phylogenyTreeNode;
        }
    }
    
    /**
     * Create a test tree
     * @return
     *          the test tree
     */
    @SuppressWarnings("all")
    protected static PhylogenyTreeNode createTestTree()
    {
        PhylogenyTreeNode phylogenyTree = new PhylogenyTreeNode();
        phylogenyTree.getStrains().add("parent");
        
        PhylogenyTreeNode child = new PhylogenyTreeNode();
        child.getStrains().add("child strain1");
        PhylogenyTreeEdge phylogenyTreeEdge = new PhylogenyTreeEdge(
                new BitSet(0),
                child,
                1.0);
        phylogenyTree.getChildEdges().add(phylogenyTreeEdge);
        
        PhylogenyTreeNode child2 = new PhylogenyTreeNode();
//        child2.getStrains().add("child strain2");
        PhylogenyTreeEdge phylogenyTreeEdge2 = new PhylogenyTreeEdge(
                new BitSet(0),
                child2,
                1.0);
        phylogenyTree.getChildEdges().add(phylogenyTreeEdge2);
        
        PhylogenyTreeNode child3 = new PhylogenyTreeNode();
        child3.getStrains().add("child strain3");
        PhylogenyTreeEdge phylogenyTreeEdge3 = new PhylogenyTreeEdge(
                new BitSet(0),
                child3,
                1.0);
        phylogenyTree.getChildEdges().add(phylogenyTreeEdge3);
        
        PhylogenyTreeNode child4 = new PhylogenyTreeNode();
        child4.getStrains().add("child strain4");
        PhylogenyTreeEdge phylogenyTreeEdge4 = new PhylogenyTreeEdge(
                new BitSet(0),
                child4,
                1.0);
        child3.getChildEdges().add(phylogenyTreeEdge4);
        
        PhylogenyTreeNode child5 = new PhylogenyTreeNode();
        child5.getStrains().add("child strain5");
        PhylogenyTreeEdge phylogenyTreeEdge5 = new PhylogenyTreeEdge(
                new BitSet(0),
                child5,
                1.0);
        child3.getChildEdges().add(phylogenyTreeEdge5);
        
        PhylogenyTreeNode child8 = null;
//        for(int i = 0; i < 12; i++)
        for(int i = 0; i < 6; i++)
        {
            child8 = new PhylogenyTreeNode();
            child8.getStrains().add("child strain6." + i);
            PhylogenyTreeEdge phylogenyTreeEdge8 = new PhylogenyTreeEdge(
                    new BitSet(0),
                    child8,
                    1.0);
            child2.getChildEdges().add(phylogenyTreeEdge8);
        }
        
//        for(int i = 0; i < 15; i++)
        for(int i = 0; i < 10; i++)
        {
            PhylogenyTreeNode child9 = new PhylogenyTreeNode();
            child9.getStrains().add("child strain7." + i);
            PhylogenyTreeEdge phylogenyTreeEdge8 = new PhylogenyTreeEdgeWithRealValue(
                    new BitSet(0),
                    child9,
                    1.0,
                    Math.random());
            child8.getChildEdges().add(phylogenyTreeEdge8);
            for(int j = 0; j < 10; j++)
            {
                child9.getStrains().add("strain8." + j);
            }
        }
        
        return phylogenyTree;
    }
    
    /**
     * A tester main
     * @param args
     *          don't care
     */
    @SuppressWarnings("all")
    public static void main(String[] args)
    {
        final PhylogenyTreeImageFactory colorFactory = new SimplePhylogenyTreeImageFactory(
                new SmoothPaintScale(
                        0.0,
                        1.0,
                        Color.BLUE,
                        Color.RED));
        JPanel colorTreePanel = new JPanel()
        {
            /**
             * {@inheritDoc}
             */
            @Override
            protected void paintComponent(Graphics g)
            {
                BufferedImage image = colorFactory.createPhylogenyImage(
                        SimplePhylogenyTreeImageFactory.createTestTree(),
                        this.getWidth(),
                        this.getHeight());
                g.drawImage(image, 0, 0, this);
            }
        };
        
        final PhylogenyTreeImageFactory noColorFactory = new SimplePhylogenyTreeImageFactory();
        JPanel noColorTreePanel = new JPanel()
        {
            /**
             * {@inheritDoc}
             */
            @Override
            protected void paintComponent(Graphics g)
            {
                BufferedImage image = noColorFactory.createPhylogenyImage(
                        SimplePhylogenyTreeImageFactory.createTestTree(),
                        this.getWidth(),
                        this.getHeight());
                g.drawImage(image, 0, 0, this);
            }
        };
        
        JFrame treeFrame = new JFrame();
        treeFrame.getContentPane().setLayout(new GridLayout());
        treeFrame.getContentPane().add(colorTreePanel);
        treeFrame.getContentPane().add(noColorTreePanel);
        treeFrame.setVisible(true);
    }
}
