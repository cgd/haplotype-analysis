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

package org.jax.haplotype.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.inference.TTest;
import org.apache.commons.math.stat.inference.TTestImpl;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeEdge;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeEdgeWithRealValue;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;
import org.jax.util.math.StatisticUtilities;

/**
 * A phylogeny significance tester
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class PhylogenySignificanceTester
{
    private static final Logger LOG = Logger.getLogger(
            PhylogenySignificanceTester.class.getName());
    
    private final TTest tTest = new TTestImpl();
    
    /**
     * Test the given phylogeny tree's edges
     * @param phylogeny
     *          the phylogeny to test
     * @param phenotypeData
     *          the strain name to phenotype data map
     * @return
     *          the significance values for the edges there will be one test
     *          value for every 
     */
    public PhylogenyTreeNode testMultipleResponseSignificance(
            PhylogenyTreeNode phylogeny,
            Map<String, List<Double>> phenotypeData)
    {
        // check that we have a phenotype for each phylogeny
        List<String> allPhyloStrains = phylogeny.getAllStrains();
        if(allPhyloStrains.size() != phenotypeData.size() ||
           !phenotypeData.keySet().containsAll(allPhyloStrains))
        {
            LOG.severe(
                    "Strain miss-match: " + allPhyloStrains.size() + " vs " +
                    phenotypeData.size());
            List<String> phenoStrains = new ArrayList<String>(
                    phenotypeData.keySet());
            Collections.sort(phenoStrains);
            Collections.sort(allPhyloStrains);
            StringBuffer phenoBuff = new StringBuffer();
            for(String strain: phenoStrains)
            {
                phenoBuff.append(strain);
                phenoBuff.append(',');
            }
            
            StringBuffer phyloBuff = new StringBuffer();
            for(String strain: allPhyloStrains)
            {
                phyloBuff.append(strain);
                phyloBuff.append(',');
            }
            
            LOG.severe("Pheno Strains: " + phenoBuff.toString());
            LOG.severe("Phylo Strains: " + phyloBuff.toString());
            LOG.severe("Phylo Tree:    " + phylogeny.toNewickFormat());
            
            throw new IllegalArgumentException(
                    "the strains in the phylogeny tree and the phenotype " +
                    "data do not match up");
        }
        
        List<PhylogenyTreeEdge> childEdgesWithPValue =
            new ArrayList<PhylogenyTreeEdge>(
                    phylogeny.getChildEdges().size());
        for(PhylogenyTreeEdge edge: phylogeny.getChildEdges())
        {
            childEdgesWithPValue.add(this.testMultipleResponseSignificanceRecursive(
                    edge,
                    phenotypeData));
        }
        PhylogenyTreeNode newNode = new PhylogenyTreeNode(
                childEdgesWithPValue,
                phylogeny.getStrains());
        return newNode;
    }

    /**
     * Recursive function that fills in the list of significance values
     * @param phylogenyEdge
     *          the edge to test
     * @return
     *          the cumulative strain set
     */
    private PhylogenyTreeEdgeWithRealValue testMultipleResponseSignificanceRecursive(
            PhylogenyTreeEdge phylogenyEdge,
            Map<String, List<Double>> phenotypeData)
    {
        PhylogenyTreeNode node = phylogenyEdge.getNode();
        
        Set<String> cumulativeBranchStrains = new HashSet<String>();
        List<PhylogenyTreeEdge> childEdgesWithPValue =
            new ArrayList<PhylogenyTreeEdge>(
                    node.getChildEdges().size());
        for(PhylogenyTreeEdge edge: node.getChildEdges())
        {
            PhylogenyTreeEdgeWithRealValue newEdgeWithPValue =
                this.testMultipleResponseSignificanceRecursive(
                        edge,
                        phenotypeData);
            childEdgesWithPValue.add(newEdgeWithPValue);
            
            // TODO make me more efficient
            // This is really inefficient. We're repeating a lot of
            // work that was done in the recursive calls here...
            cumulativeBranchStrains.addAll(
                    newEdgeWithPValue.getNode().getAllStrains());
        }
        
        cumulativeBranchStrains.addAll(node.getStrains());
        
        double significanceValue = 1.0;
        int inPhyloCount = cumulativeBranchStrains.size();
        if(inPhyloCount >= 2 && (phenotypeData.size() - inPhyloCount) >= 2)
        {
            DescriptiveStatistics inPhyloStats = new DescriptiveStatistics();
            DescriptiveStatistics outPhyloStats = new DescriptiveStatistics();
            
            for(Entry<String, List<Double>> entry: phenotypeData.entrySet())
            {
                if(cumulativeBranchStrains.contains(entry.getKey()))
                {
                    inPhyloStats.addValue(StatisticUtilities.calculateMean(
                            entry.getValue()));
                }
                else
                {
                    outPhyloStats.addValue(StatisticUtilities.calculateMean(
                            entry.getValue()));
                }
            }
            
            try
            {
                significanceValue = this.tTest.tTest(inPhyloStats, outPhyloStats);
            }
            catch(Exception ex)
            {
                LOG.log(Level.SEVERE,
                        "t test failed",
                        ex);
            }
        }
        PhylogenyTreeNode newNode = new PhylogenyTreeNode(
                childEdgesWithPValue,
                node.getStrains());
        return new PhylogenyTreeEdgeWithRealValue(
                    phylogenyEdge.getSdpBits(),
                    newNode,
                    phylogenyEdge.getEdgeLength(),
                    significanceValue);
    }
}
