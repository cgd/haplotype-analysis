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

package org.jax.haplotype.analysis;

/**
 * For performing an association test using native EMMA library
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class EMMAAssociationTest
{
    static
    {
        System.loadLibrary("emma");
    }
    
    /**
     * Native function for performing an EMMA scan
     * @param strainCount   the number of strains
     * @param phenos        the (flattened) phenotype matrix
     * @param genos         the (flattened) genotype matrix
     * @param kinship       the (flattened) kinship matrix. this can be null
     *                      in which case it will be estimated from the given
     *                      genotypes
     * @return              the resulting (flattened) pvalue matrix
     */
    private static native double[] emmaScan(
            int strainCount,
            double[] phenos,
            double[] genos,
            double[] kinship);
    
    /**
     * Native function for performing an EMMA scan
     * @param strainCount   the number of strains
     * @param genos         the (flattened) genotype matrix
     * @return              the resulting (flattened) kinship matrix
     */
    private static native double[] calculateKinship(
            int strainCount,
            double[] genos);
}
