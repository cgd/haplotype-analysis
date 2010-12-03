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

package org.jax.haplotype.analysis.experimentdesign;

import java.io.BufferedReader;
import java.io.FileReader;

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
    
    private static class INeedThisLameClassBecauseJavaHasNoTupleSupport
    {
        public double[] values;
        public int columnCount;
    }
    
    private static INeedThisLameClassBecauseJavaHasNoTupleSupport readDoubleMatrix(String fileName) throws Exception
    {
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        String currLine = null;
        int columnCount = -1;
        int valCount = 0;
        while((currLine = reader.readLine()) != null)
        {
            String[] valStrings = currLine.split("\\s+");
            if(columnCount == -1)
            {
                columnCount = valStrings.length;
            }
            else if(valStrings.length != columnCount)
            {
                throw new Exception("Count missmatch!");
            }
            
            for(String valStr: valStrings)
            {
                valCount++;
            }
        }
        
        reader = new BufferedReader(new FileReader(fileName));
        double[] doubles = new double[valCount];
        valCount = 0;
        while((currLine = reader.readLine()) != null)
        {
            String[] valStrings = currLine.split("\\s+");
            for(String valStr: valStrings)
            {
                if(valStr.equals("NA"))
                {
                    doubles[valCount] = Double.NaN;
                }
                else
                {
                    doubles[valCount] = Double.parseDouble(valStr);
                }
                
                valCount++;
            }
        }
        
        INeedThisLameClassBecauseJavaHasNoTupleSupport lame =
            new INeedThisLameClassBecauseJavaHasNoTupleSupport();
        lame.values = doubles;
        lame.columnCount = columnCount;
        
        return lame;
    }
    
//    /**
//     * Main function for testing
//     * @param args  don't care
//     * @throws Exception
//     */
//    public static void main(String[] args) throws Exception
//    {
////        java.util.Random rand = new java.util.Random();
////        
////        int strainCount = 20;
////        int phenoCount = 2;
////        int snpCount = 50;
////        
////        double[] phenos = new double[strainCount * phenoCount];
////        for(int i = 0; i < phenos.length; i++)
////        {
////            phenos[i] = rand.nextDouble();
////        }
////        
////        double[] genos = new double[snpCount * strainCount];
////        for(int i = 0; i < genos.length; i++)
////        {
////            genos[i] = rand.nextDouble();
////        }
////        
////        double[] result = EMMAAssociationTest.emmaScan(strainCount, phenos, genos);
////        System.out.println(SequenceUtilities.toString(SequenceUtilities.toDoubleList(result)));
//        
////        INeedThisLameClassBecauseJavaHasNoTupleSupport hdlLame = readDoubleMatrix("data/hmdp.hdl.tab");
////        INeedThisLameClassBecauseJavaHasNoTupleSupport snpLame = readDoubleMatrix("data/hmdp_1000.snps");
////        INeedThisLameClassBecauseJavaHasNoTupleSupport kinshipLame = readDoubleMatrix("data/hmdp.K");
////        
////        if(hdlLame.columnCount != snpLame.columnCount)
////        {
////            throw new Exception("Column counts don't match!");
////        }
////        
////        if(kinshipLame.values.length != (snpLame.columnCount * snpLame.columnCount))
////        {
////            throw new Exception("kinship size doesn't match: " + kinshipLame.values.length + " vs " + (snpLame.columnCount * snpLame.columnCount));
////        }
////        
////        double[] scanResult = emmaScan(hdlLame.columnCount, hdlLame.values, snpLame.values, kinshipLame.values);
////        
////        for(double pValue: scanResult)
////        {
////            System.out.println(pValue);
////        }
//        
//        INeedThisLameClassBecauseJavaHasNoTupleSupport snpLame = readDoubleMatrix("data/hmdp.snps");
//        double[] kinship = calculateKinship(snpLame.columnCount, snpLame.values);
//        for(int i = 0; i < kinship.length; i++)
//        {
//            System.out.print(kinship[i]);
//            
//            if((i + 1) % snpLame.columnCount == 0)
//            {
//                System.out.println();
//            }
//            else
//            {
//                System.out.print('\t');
//            }
//        }
//    }
}
