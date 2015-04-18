/**
 * IBM Confidential
 * OCO Source Materials
 * (C) Copyright IBM Corp. 2010, 2015
 * The source code for this program is not published or otherwise divested of its trade secrets, irrespective of what has been deposited with the U.S. Copyright Office.
 */

package com.ibm.bi.dml.runtime.instructions.mr;

import java.util.ArrayList;

import com.ibm.bi.dml.lops.MapMultChain.ChainType;
import com.ibm.bi.dml.runtime.DMLRuntimeException;
import com.ibm.bi.dml.runtime.DMLUnsupportedOperationException;
import com.ibm.bi.dml.runtime.instructions.Instruction;
import com.ibm.bi.dml.runtime.instructions.InstructionUtils;
import com.ibm.bi.dml.runtime.matrix.data.MatrixBlock;
import com.ibm.bi.dml.runtime.matrix.data.MatrixIndexes;
import com.ibm.bi.dml.runtime.matrix.data.MatrixValue;
import com.ibm.bi.dml.runtime.matrix.mapred.CachedValueMap;
import com.ibm.bi.dml.runtime.matrix.mapred.DistributedCacheInput;
import com.ibm.bi.dml.runtime.matrix.mapred.IndexedMatrixValue;
import com.ibm.bi.dml.runtime.matrix.mapred.MRBaseForCommonInstructions;

/**
 * 
 */
public class MapMultChainInstruction extends MRInstruction 
{
	@SuppressWarnings("unused")
	private static final String _COPYRIGHT = "Licensed Materials - Property of IBM\n(C) Copyright IBM Corp. 2010, 2015\n" +
                                             "US Government Users Restricted Rights - Use, duplication  disclosure restricted by GSA ADP Schedule Contract with IBM Corp.";
	
	private ChainType _chainType = null;
	
	private byte _input1 = -1;
	private byte _input2 = -1;
	private byte _input3 = -1;	
	
	/**
	 * Two matrix inputs - type XtXv
	 * 
	 * @param type
	 * @param in1
	 * @param in2
	 * @param out
	 * @param istr
	 */
	public MapMultChainInstruction(ChainType type, byte in1, byte in2, byte out, String istr)
	{
		super(null, out);
		
		_chainType = type;
		
		_input1 = in1;
		_input2 = in2;
		_input3 = -1;
		
		mrtype = MRINSTRUCTION_TYPE.MapMultChain;
		instString = istr;
	}

	/**
	 * Three matrix inputs - type XtwXv
	 * 
	 * @param type
	 * @param in1
	 * @param in2
	 * @param in3
	 * @param out
	 * @param istr
	 */
	public MapMultChainInstruction(ChainType type, byte in1, byte in2, byte in3, byte out, String istr)
	{
		super(null, out);
		
		_chainType = type;
		
		_input1 = in1;
		_input2 = in2;
		_input3 = in3;
		
		mrtype = MRINSTRUCTION_TYPE.MapMultChain;
		instString = istr;
	}
	
	public ChainType getChainType()
	{
		return _chainType;
	}
	
	public byte getInput1() {
		return _input1;
	}

	public byte getInput2() {
		return _input2;
	}

	public byte getInput3() {
		return _input3;
	}

	/**
	 * 
	 * @param str
	 * @return
	 * @throws DMLRuntimeException
	 */
	public static Instruction parseInstruction( String str ) 
		throws DMLRuntimeException 
	{		
		//check number of fields (2/3 inputs, output, type)
		InstructionUtils.checkNumFields ( str, 4, 5 );
		
		//parse instruction parts (without exec type)
		String[] parts = InstructionUtils.getInstructionParts( str );		
		byte in1 = Byte.parseByte(parts[1]);
		byte in2 = Byte.parseByte(parts[2]);
		
		if( parts.length==5 )
		{
			byte out = Byte.parseByte(parts[3]);
			ChainType type = ChainType.valueOf(parts[4]);
			
			return new MapMultChainInstruction(type, in1, in2, out, str);
		}
		else //parts.length==6
		{
			byte in3 = Byte.parseByte(parts[3]);
			byte out = Byte.parseByte(parts[4]);
			ChainType type = ChainType.valueOf(parts[5]);
		
			return new MapMultChainInstruction(type, in1, in2, in3, out, str);
		}	
	}
	
	/**
	 * Determines if the given index is only used via distributed cache in
	 * the given instruction string (used during setup of distributed cache
	 * to detect redundant job inputs).
	 * 
	 * @param inst
	 * @param index
	 * @return
	 */
	public static boolean isDistCacheOnlyIndex( String inst, byte index )
	{
		boolean ret = false;
		
		//parse instruction parts (with exec type)
		String[] parts = inst.split(Instruction.OPERAND_DELIM);
		if( parts.length==6 ){
			byte in1 = Byte.parseByte(parts[2].split(Instruction.DATATYPE_PREFIX)[0]);
			byte in2 = Byte.parseByte(parts[3].split(Instruction.DATATYPE_PREFIX)[0]);
			ret = (index==in2 && index!=in1);
		}
		else if( parts.length==7 ){
			byte in1 = Byte.parseByte(parts[2].split(Instruction.DATATYPE_PREFIX)[0]);
			byte in2 = Byte.parseByte(parts[3].split(Instruction.DATATYPE_PREFIX)[0]);
			byte in3 = Byte.parseByte(parts[4].split(Instruction.DATATYPE_PREFIX)[0]);
			ret = (index==in2 && index!=in1) || (index==in3 && index!=in1);
		}
		
		return ret;
	}
	
	public static void addDistCacheIndex( String inst, ArrayList<Byte> indexes )
	{
		//parse instruction parts (with exec type)
		String[] parts = inst.split(Instruction.OPERAND_DELIM);
		if( parts.length==6 ){
			byte in2 = Byte.parseByte(parts[3].split(Instruction.DATATYPE_PREFIX)[0]);
			indexes.add(in2);
		}
		else if( parts.length==7 ){
			byte in2 = Byte.parseByte(parts[3].split(Instruction.DATATYPE_PREFIX)[0]);
			byte in3 = Byte.parseByte(parts[4].split(Instruction.DATATYPE_PREFIX)[0]);
			indexes.add(in2);
			indexes.add(in3);
		}
	}
	
	@Override
	public byte[] getInputIndexes() 
	{
		if( _chainType==ChainType.XtXv )
			return new byte[]{_input1, _input2};
		else
			return new byte[]{_input1, _input2, _input3};
	}

	@Override
	public byte[] getAllIndexes() 
	{
		if( _chainType==ChainType.XtXv )
			return new byte[]{_input1, _input2, output};
		else
			return new byte[]{_input1, _input2, _input3, output};
	}
	

	@Override
	public void processInstruction(Class<? extends MatrixValue> valueClass, CachedValueMap cachedValues, 
			           IndexedMatrixValue tempValue, IndexedMatrixValue zeroInput, int blockRowFactor, int blockColFactor)
		throws DMLUnsupportedOperationException, DMLRuntimeException 
	{
		ArrayList<IndexedMatrixValue> blkList = cachedValues.get(_input1);
		if( blkList !=null )
			for(IndexedMatrixValue imv : blkList)
			{
				if(imv==null)
					continue;
				MatrixIndexes inIx = imv.getIndexes();
				MatrixValue inVal = imv.getValue();
				
				//allocate space for the output value
				IndexedMatrixValue iout = null;
				if(output==_input1)
					iout=tempValue;
				else
					iout=cachedValues.holdPlace(output, valueClass);
				
				MatrixIndexes outIx = iout.getIndexes();
				MatrixValue outVal = iout.getValue();
				
				//process instruction
				if( _chainType == ChainType.XtXv )
					processXtXvOperations(inIx, inVal, outIx, outVal);
				else
					processXtwXvOperations(inIx, inVal, outIx, outVal);
				
				//put the output value in the cache
				if(iout==tempValue)
					cachedValues.add(output, iout);
			}
	}

	
	/**
	 * Chain implementation for r = (t(X)%*%(X%*%v))
	 * (implemented as r = (t(t(X%*%v)%*%X))
	 * 
	 * @param inIx
	 * @param inVal
	 * @param outIx
	 * @param outVal
	 * @throws DMLRuntimeException 
	 * @throws DMLUnsupportedOperationException 
	 */
	private void processXtXvOperations(MatrixIndexes inIx, MatrixValue inVal, MatrixIndexes outIx, MatrixValue outVal ) 
		throws DMLRuntimeException, DMLUnsupportedOperationException
	{
		DistributedCacheInput dcInput2 = MRBaseForCommonInstructions.dcValues.get(_input2); //v
		MatrixBlock Xi = (MatrixBlock)inVal;
		MatrixBlock v = (MatrixBlock) dcInput2.getDataBlock(1, 1).getValue();
		
		//process core block operation
		Xi.chainMatrixMultOperations(v, null, (MatrixBlock) outVal, ChainType.XtXv);
		outIx.setIndexes(1, 1);
	}
	
	/**
	 * Chain implementation for r = (t(X)%*%(w*(X%*%v)))
	 * (implemented as r = (t(t((X%*%v)*w)%*%X))
	 * 
	 * @param inIx
	 * @param inVal
	 * @param outIx
	 * @param outVal
	 * @throws DMLRuntimeException 
	 * @throws DMLUnsupportedOperationException 
	 */
	private void processXtwXvOperations(MatrixIndexes inIx, MatrixValue inVal, MatrixIndexes outIx, MatrixValue outVal )
		throws DMLRuntimeException, DMLUnsupportedOperationException
	{
		DistributedCacheInput dcInput2 = MRBaseForCommonInstructions.dcValues.get(_input2); //v
		DistributedCacheInput dcInput3 = MRBaseForCommonInstructions.dcValues.get(_input3); //w
		MatrixBlock Xi = (MatrixBlock) inVal;
		MatrixBlock v = (MatrixBlock) dcInput2.getDataBlock(1, 1).getValue();
		MatrixBlock w = (MatrixBlock) dcInput3.getDataBlock((int)inIx.getRowIndex(), 1).getValue();
		
		//process core block operation
		Xi.chainMatrixMultOperations(v, w, (MatrixBlock) outVal, ChainType.XtwXv);
		outIx.setIndexes(1, 1);
	}
}
