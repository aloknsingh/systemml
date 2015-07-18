/**
 * IBM Confidential
 * OCO Source Materials
 * (C) Copyright IBM Corp. 2010, 2015
 * The source code for this program is not published or otherwise divested of its trade secrets, irrespective of what has been deposited with the U.S. Copyright Office.
 */

package com.ibm.bi.dml.test.integration.functions.tertiary;

import java.util.HashMap;

import org.junit.Test;

import com.ibm.bi.dml.api.DMLScript;
import com.ibm.bi.dml.api.DMLScript.RUNTIME_PLATFORM;
import com.ibm.bi.dml.lops.LopProperties.ExecType;
import com.ibm.bi.dml.runtime.matrix.data.MatrixValue.CellIndex;
import com.ibm.bi.dml.test.integration.AutomatedTestBase;
import com.ibm.bi.dml.test.integration.TestConfiguration;
import com.ibm.bi.dml.test.utils.TestUtils;

/**
 * 
 */
public class CentralMomentWeightsTest extends AutomatedTestBase 
{
	@SuppressWarnings("unused")
	private static final String _COPYRIGHT = "Licensed Materials - Property of IBM\n(C) Copyright IBM Corp. 2010, 2015\n" +
                                             "US Government Users Restricted Rights - Use, duplication  disclosure restricted by GSA ADP Schedule Contract with IBM Corp.";
	
	private final static String TEST_NAME = "CentralMomentWeights";
	private final static String TEST_DIR = "functions/tertiary/";
	private final static double eps = 1e-10;
	
	private final static int rows = 1871;
	private final static int maxVal = 7; 
	private final static double sparsity1 = 0.65;
	private final static double sparsity2 = 0.05;
	
	
	@Override
	public void setUp() 
	{
		TestUtils.clearAssertionInformation();
		addTestConfiguration(TEST_NAME, new TestConfiguration(TEST_DIR, TEST_NAME, new String[] { "R" })   ); 
	}

	@Test
	public void testCentralMoment2WeightsDenseCP() 
	{
		runCentralMomentTest(2, false, ExecType.CP);
	}
	
	@Test
	public void testCentralMoment3WeightsDenseCP() 
	{
		runCentralMomentTest(3, false, ExecType.CP);
	}
	
	@Test
	public void testCentralMoment4WeightsDenseCP() 
	{
		runCentralMomentTest(4, false, ExecType.CP);
	}
	
	@Test
	public void testCentralMoment2WeightsSparseCP() 
	{
		runCentralMomentTest(2, true, ExecType.CP);
	}
	
	@Test
	public void testCentralMoment3WeightsSparseCP() 
	{
		runCentralMomentTest(3, true, ExecType.CP);
	}
	
	@Test
	public void testCentralMoment4WeightsSparseCP() 
	{
		runCentralMomentTest(4, true, ExecType.CP);
	}
	
	@Test
	public void testCentralMoment2WeightsDenseMR() 
	{
		runCentralMomentTest(2, false, ExecType.MR);
	}
	
	@Test
	public void testCentralMoment3WeightsDenseMR() 
	{
		runCentralMomentTest(3, false, ExecType.MR);
	}
	
	@Test
	public void testCentralMoment4WeightsDenseMR() 
	{
		runCentralMomentTest(4, false, ExecType.MR);
	}
	
	@Test
	public void testCentralMoment2WeightsSparseMR() 
	{
		runCentralMomentTest(2, true, ExecType.MR);
	}
	
	@Test
	public void testCentralMoment3WeightsSparseMR() 
	{
		runCentralMomentTest(3, true, ExecType.MR);
	}
	
	@Test
	public void testCentralMoment4WeightsSparseMR() 
	{
		runCentralMomentTest(4, true, ExecType.MR);
	}
	
	@Test
	public void testCentralMoment2WeightsDenseSP() 
	{
		runCentralMomentTest(2, false, ExecType.SPARK);
	}
	
	@Test
	public void testCentralMoment3WeightsDenseSP() 
	{
		runCentralMomentTest(3, false, ExecType.SPARK);
	}
	
	@Test
	public void testCentralMoment4WeightsDenseSP() 
	{
		runCentralMomentTest(4, false, ExecType.SPARK);
	}
	
	@Test
	public void testCentralMoment2WeightsSparseSP() 
	{
		runCentralMomentTest(2, true, ExecType.SPARK);
	}
	
	@Test
	public void testCentralMoment3WeightsSparseSP() 
	{
		runCentralMomentTest(3, true, ExecType.SPARK);
	}
	
	@Test
	public void testCentralMoment4WeightsSparseSP() 
	{
		runCentralMomentTest(4, true, ExecType.SPARK);
	}
	
	/**
	 * 
	 * @param sparseM1
	 * @param sparseM2
	 * @param instType
	 */
	private void runCentralMomentTest( int order, boolean sparse, ExecType et)
	{
		//rtplatform for MR
		RUNTIME_PLATFORM platformOld = rtplatform;
		switch( et ){
			case MR: rtplatform = RUNTIME_PLATFORM.HADOOP; break;
			case SPARK: rtplatform = RUNTIME_PLATFORM.SPARK; break;
			default: rtplatform = RUNTIME_PLATFORM.HYBRID; break;
		}
	
		boolean sparkConfigOld = DMLScript.USE_LOCAL_SPARK_CONFIG;
		if( rtplatform == RUNTIME_PLATFORM.SPARK )
			DMLScript.USE_LOCAL_SPARK_CONFIG = true;
		
		try
		{
			TestConfiguration config = getTestConfiguration(TEST_NAME);
			
			String HOME = SCRIPT_DIR + TEST_DIR;
			fullDMLScriptName = HOME + TEST_NAME + ".dml";
			programArgs = new String[]{"-args", HOME + INPUT_DIR + "A",
					                            HOME + INPUT_DIR + "B",
					                        Integer.toString(order),
					                        HOME + OUTPUT_DIR + "R"};
			fullRScriptName = HOME + TEST_NAME + ".R";
			rCmd = "Rscript" + " " + fullRScriptName + " " + 
			       HOME + INPUT_DIR + " " + order + " "+ HOME + EXPECTED_DIR;
			
			loadTestConfiguration(config);
	
			//generate actual dataset (always dense because values <=0 invalid)
			double sparsitya = sparse ? sparsity2 : sparsity1;
			double[][] A = getRandomMatrix(rows, 1, 1, maxVal, sparsitya, 7); 
			writeInputMatrixWithMTD("A", A, true);
			
			double[][] B = getRandomMatrix(rows, 1, 1, 1, 1.0, 34); 
			writeInputMatrixWithMTD("B", B, true);	
			
			runTest(true, false, null, -1); 
			runRScript(true); 
			
			//compare matrices 
			HashMap<CellIndex, Double> dmlfile = readDMLMatrixFromHDFS("R");
			HashMap<CellIndex, Double> rfile  = readRMatrixFromFS("R");
			TestUtils.compareMatrices(dmlfile, rfile, eps, "Stat-DML", "Stat-R");
		}
		finally
		{
			rtplatform = platformOld;
			DMLScript.USE_LOCAL_SPARK_CONFIG = sparkConfigOld;
		}
	}

}