/**
 * (C) Copyright IBM Corp. 2010, 2015
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

package com.ibm.bi.dml.api;

import java.util.ArrayList;
import java.util.HashMap;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.mllib.linalg.DenseVector;
import org.apache.spark.sql.DataFrame;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.RowFactory;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.types.StructType;

import scala.Tuple2;

import com.ibm.bi.dml.parser.DMLTranslator;
import com.ibm.bi.dml.runtime.DMLRuntimeException;
import com.ibm.bi.dml.runtime.instructions.spark.functions.ConvertMatrixBlockToIJVLines;
import com.ibm.bi.dml.runtime.instructions.spark.functions.GetMLBlock;
import com.ibm.bi.dml.runtime.instructions.spark.utils.RDDConverterUtilsExt;
import com.ibm.bi.dml.runtime.matrix.MatrixCharacteristics;
import com.ibm.bi.dml.runtime.matrix.data.MatrixBlock;
import com.ibm.bi.dml.runtime.matrix.data.MatrixIndexes;
import com.ibm.bi.dml.runtime.util.UtilFunctions;

/**
 * This is a simple container object that returns the output of execute from MLContext 
 *
 */
public class MLOutput {
	
	
	
	HashMap<String, JavaPairRDD<MatrixIndexes,MatrixBlock>> _outputs;
	private HashMap<String, MatrixCharacteristics> _outMetadata = null;
	
	public MLOutput(HashMap<String, JavaPairRDD<MatrixIndexes,MatrixBlock>> outputs, HashMap<String, MatrixCharacteristics> outMetadata) {
		this._outputs = outputs;
		this._outMetadata = outMetadata;
	}
	
	public JavaPairRDD<MatrixIndexes,MatrixBlock> getBinaryBlockedRDD(String varName) throws DMLRuntimeException {
		if(_outputs.containsKey(varName)) {
			return _outputs.get(varName);
		}
		throw new DMLRuntimeException("Variable " + varName + " not found in the output symbol table.");
	}
	
	public MatrixCharacteristics getMatrixCharacteristics(String varName) throws DMLRuntimeException {
		if(_outputs.containsKey(varName)) {
			return _outMetadata.get(varName);
		}
		throw new DMLRuntimeException("Variable " + varName + " not found in the output symbol table.");
	}
	
	/**
	 * Note, the output DataFrame has an additional column ID.
	 * An easy way to get DataFrame without ID is by df.sort("ID").drop("ID")
	 * @param sqlContext
	 * @param varName
	 * @return
	 * @throws DMLRuntimeException
	 */
	public DataFrame getDF(SQLContext sqlContext, String varName) throws DMLRuntimeException {
		JavaPairRDD<MatrixIndexes,MatrixBlock> rdd = getBinaryBlockedRDD(varName);
		if(rdd != null) {
			MatrixCharacteristics mc = _outMetadata.get(varName);
			return RDDConverterUtilsExt.binaryBlockToDataFrame(rdd, mc, sqlContext);
		}
		throw new DMLRuntimeException("Variable " + varName + " not found in the output symbol table.");
	}
	
	/**
	 * 
	 * @param sqlContext
	 * @param varName
	 * @param outputVector if true, returns DataFrame with two column: ID and org.apache.spark.mllib.linalg.Vector
	 * @return
	 * @throws DMLRuntimeException
	 */
	public DataFrame getDF(SQLContext sqlContext, String varName, boolean outputVector) throws DMLRuntimeException {
		if(outputVector) {
			JavaPairRDD<MatrixIndexes,MatrixBlock> rdd = getBinaryBlockedRDD(varName);
			if(rdd != null) {
				MatrixCharacteristics mc = _outMetadata.get(varName);
				return RDDConverterUtilsExt.binaryBlockToVectorDataFrame(rdd, mc, sqlContext);
			}
			throw new DMLRuntimeException("Variable " + varName + " not found in the output symbol table.");
		}
		else {
			return getDF(sqlContext, varName);
		}
		
	}
	
	public JavaRDD<String> getStringRDD(String varName, String format) throws DMLRuntimeException {
		if(format.compareTo("text") == 0) {
			JavaPairRDD<MatrixIndexes,MatrixBlock> in1 = getBinaryBlockedRDD(varName);
			JavaRDD<String> ijv = in1.flatMap(new ConvertMatrixBlockToIJVLines(DMLTranslator.DMLBlockSize, DMLTranslator.DMLBlockSize));
			return ijv;
		}
//		else if(format.compareTo("csv") == 0) {
//			
//		}
		else {
			throw new DMLRuntimeException("The output format:" + format + " is not implemented yet.");
		}
		
	}
	
	public MLMatrix getMLMatrix(MLContext ml, SQLContext sqlContext, String varName) throws DMLRuntimeException {
		JavaPairRDD<MatrixIndexes,MatrixBlock> rdd = getBinaryBlockedRDD(varName);
		if(rdd != null) {
			MatrixCharacteristics mc = getMatrixCharacteristics(varName);
			StructType schema = MLBlock.getDefaultSchemaForBinaryBlock();
			return new MLMatrix(sqlContext.createDataFrame(rdd.map(new GetMLBlock()).rdd(), schema), mc, ml);
		}
		throw new DMLRuntimeException("Variable " + varName + " not found in the output symbol table.");
	}
	
//	/**
//	 * Experimental: Please use this with caution as it will fail in many corner cases.
//	 * @return org.apache.spark.mllib.linalg.distributed.BlockMatrix
//	 * @throws DMLRuntimeException 
//	 */
//	public BlockMatrix getMLLibBlockedMatrix(MLContext ml, SQLContext sqlContext, String varName) throws DMLRuntimeException {
//		return getMLMatrix(ml, sqlContext, varName).toBlockedMatrix();
//	}
	
	public static class ProjectRows implements PairFlatMapFunction<Tuple2<MatrixIndexes,MatrixBlock>, Long, Tuple2<Long, Double[]>> {
		private static final long serialVersionUID = -4792573268900472749L;
		long rlen; long clen;
		int brlen; int bclen;
		public ProjectRows(long rlen, long clen, int brlen, int bclen) {
			this.rlen = rlen;
			this.clen = clen;
			this.brlen = brlen;
			this.bclen = bclen;
		}

		@Override
		public Iterable<Tuple2<Long, Tuple2<Long, Double[]>>> call(Tuple2<MatrixIndexes, MatrixBlock> kv) throws Exception {
			// ------------------------------------------------------------------
    		//	Compute local block size: 
    		// Example: For matrix: 1500 X 1100 with block length 1000 X 1000
    		// We will have four local block sizes (1000X1000, 1000X100, 500X1000 and 500X1000)
    		long blockRowIndex = kv._1.getRowIndex();
    		long blockColIndex = kv._1.getColumnIndex();
    		int lrlen = UtilFunctions.computeBlockSize(rlen, blockRowIndex, brlen);
    		int lclen = UtilFunctions.computeBlockSize(clen, blockColIndex, bclen);
    		// ------------------------------------------------------------------
			
			long startRowIndex = (kv._1.getRowIndex()-1) * bclen;
			MatrixBlock blk = kv._2;
			ArrayList<Tuple2<Long, Tuple2<Long, Double[]>>> retVal = new ArrayList<Tuple2<Long,Tuple2<Long,Double[]>>>();
			for(int i = 0; i < lrlen; i++) {
				Double[] partialRow = new Double[lclen];
				for(int j = 0; j < lclen; j++) {
					partialRow[j] = blk.getValue(i, j);
				}
				retVal.add(new Tuple2<Long, Tuple2<Long,Double[]>>(startRowIndex + i, new Tuple2<Long,Double[]>(kv._1.getColumnIndex(), partialRow)));
			}
			return (Iterable<Tuple2<Long, Tuple2<Long, Double[]>>>) retVal;
		}
	}
	
	public static class ConvertDoubleArrayToRows implements Function<Tuple2<Long, Iterable<Tuple2<Long, Double[]>>>, Row> {
		private static final long serialVersionUID = 4441184411670316972L;
		
		int bclen; long clen;
		boolean outputVector;
		public ConvertDoubleArrayToRows(long clen, int bclen, boolean outputVector) {
			this.bclen = bclen;
			this.clen = clen;
			this.outputVector = outputVector;
		}

		@Override
		public Row call(Tuple2<Long, Iterable<Tuple2<Long, Double[]>>> arg0)
				throws Exception {
			
			HashMap<Long, Double[]> partialRows = new HashMap<Long, Double[]>();
			int sizeOfPartialRows = 0;
			for(Tuple2<Long, Double[]> kv : arg0._2) {
				partialRows.put(kv._1, kv._2);
				sizeOfPartialRows += kv._2.length;
			}
			
			// Insert first row as row index
			Object[] row = null;
			if(outputVector) {
				row = new Object[2];
				double [] vecVals = new double[sizeOfPartialRows];
				
				for(long columnBlockIndex = 1; columnBlockIndex <= partialRows.size(); columnBlockIndex++) {
					if(partialRows.containsKey(columnBlockIndex)) {
						Double [] array = partialRows.get(columnBlockIndex);
						// ------------------------------------------------------------------
						//	Compute local block size: 
						int lclen = UtilFunctions.computeBlockSize(clen, columnBlockIndex, bclen);
						// ------------------------------------------------------------------
						if(array.length != lclen) {
							throw new Exception("Incorrect double array provided by ProjectRows");
						}
						for(int i = 0; i < lclen; i++) {
							vecVals[(int) ((columnBlockIndex-1)*bclen + i)] = array[i];
						}
					}
					else {
						throw new Exception("The block for column index " + columnBlockIndex + " is missing. Make sure the last instruction is not returning empty blocks");
					}
				}
				
				long rowIndex = arg0._1;
				row[0] = new Double(rowIndex);
				row[1] = new DenseVector(vecVals); // breeze.util.JavaArrayOps.arrayDToDv(vecVals);
			}
			else {
				row = new Double[sizeOfPartialRows + 1];
				long rowIndex = arg0._1;
				row[0] = new Double(rowIndex);
				for(long columnBlockIndex = 1; columnBlockIndex <= partialRows.size(); columnBlockIndex++) {
					if(partialRows.containsKey(columnBlockIndex)) {
						Double [] array = partialRows.get(columnBlockIndex);
						// ------------------------------------------------------------------
						//	Compute local block size: 
						int lclen = UtilFunctions.computeBlockSize(clen, columnBlockIndex, bclen);
						// ------------------------------------------------------------------
						if(array.length != lclen) {
							throw new Exception("Incorrect double array provided by ProjectRows");
						}
						for(int i = 0; i < lclen; i++) {
							row[(int) ((columnBlockIndex-1)*bclen + i) + 1] = array[i];
						}
					}
					else {
						throw new Exception("The block for column index " + columnBlockIndex + " is missing. Make sure the last instruction is not returning empty blocks");
					}
				}
			}
			Object[] row_fields = row;
			return RowFactory.create(row_fields);
		}
		
	}
	
}
