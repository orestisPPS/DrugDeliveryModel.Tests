using System;
using System.Collections.Generic;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.DrugDeliveryModel.Tests.Commons;

public class BoundaryConditionsUtility
{
	public enum BoundaryConditionCase
	{
		LeftDirichlet,
		LeftPointFlux,
		LeftDistributedForce,
		
		RightDirichlet,
		RightPointFlux,
		RightDistributedFlux,
		
		FrontDirichlet,
		FrontPointFlux,
		FrontDistributedFlux,
		
		BackDirichlet,
		BackPointFlux,
		BackDistributedFlux,
		
		BottomDirichlet,
		BottomPointFlux,
		BottomDistributedFlux,
		
		TopDirichlet,
		TopPointFlux,
		TopDistributedFlux
	}

	public static void AssignStructuralDirichletBCsToModel(Model model,
		List<(BoundaryConditionCase, StructuralDof[], double[][], double[])> bcs,
		double tolerance)
	{
		var constraints = new List<INodalDisplacementBoundaryCondition>();
		var emptyloads = new List<INodalLoadBoundaryCondition>();
		
		foreach (var bcSet in bcs)
		{
			var bcFace = bcSet.Item1;
			var bcDofs = bcSet.Item2;
			var faceCoords = bcSet.Item3;
			var bcValues = bcSet.Item4;
			
			foreach (var node in model.NodesDictionary)
			{
				var nodalCoords = new double[] { node.Value.X, node.Value.Y, node.Value.Z };
				
				if (bcFace == BoundaryConditionCase.RightDirichlet || bcFace == BoundaryConditionCase.LeftDirichlet
				    && Math.Abs(nodalCoords[0] - faceCoords[0][0]) < tolerance)
				{
					ApplyStructuralDirichletBCToNode(constraints, node.Value, bcDofs, bcValues);
				}
				
				if (bcFace == BoundaryConditionCase.FrontDirichlet || bcFace == BoundaryConditionCase.BackDirichlet &&
				         Math.Abs(nodalCoords[1] - faceCoords[0][1]) < tolerance)
				{
					ApplyStructuralDirichletBCToNode(constraints, node.Value, bcDofs, bcValues);
				}
				
				if (bcFace == BoundaryConditionCase.TopDirichlet || bcFace == BoundaryConditionCase.BottomDirichlet &&
				         Math.Abs(nodalCoords[2] - faceCoords[0][2]) < tolerance)
				{
					ApplyStructuralDirichletBCToNode(constraints, node.Value, bcDofs, bcValues);
				}
			}
		}
		model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, emptyloads));
	}

	public static void AssignStructuralNeumannBCsToModel(Model model, List<(BoundaryConditionCase, StructuralDof[], double[][], double[])> bcs,
	    double tolerance)
    {
        var loads = new List<INodalLoadBoundaryCondition>();
        
        foreach (var node in model.NodesDictionary)
        {
        	var nodalCoords = new double[] { node.Value.X, node.Value.Y, node.Value.Z };

            foreach (var bcSet in bcs)
            {
                var loadCase = bcSet.Item1;
                var bcDofs = bcSet.Item2;
                var loadCoords = bcSet.Item3;
                var bcValues = bcSet.Item4;
				
                for (int i = 0; i <= loadCoords.GetLength(0) - 1; i++)
                {
	                var coordDistanceCondition = Math.Abs(nodalCoords[0] - loadCoords[i][0]) < tolerance &&
													 Math.Abs(nodalCoords[1] - loadCoords[i][1]) < tolerance &&
													 Math.Abs(nodalCoords[2] - loadCoords[i][2]) < tolerance;
	                
	                if (coordDistanceCondition is true)
		                ApplyStructuralNeumannBCToNode(loads, node.Value, bcDofs, bcValues);
                }
            }
        }
        var modelNeumannConditions = new StructuralBoundaryConditionSet(new List<INodalDisplacementBoundaryCondition>(), loads);
        model.BoundaryConditions.Add(modelNeumannConditions);
    }
    
	private static void ApplyStructuralNeumannBCToNode(List<INodalLoadBoundaryCondition> loads,
	    INode node, StructuralDof[] dofs, double[] values)
    {
	    for (int i = 0; i < dofs.Length; i++)
	    {
		    loads.Add(new NodalLoad(node, dofs[i], values[i]));
	    }
    }

    private static void ApplyStructuralDirichletBCToNode(List<INodalDisplacementBoundaryCondition> constraints,
		INode node, StructuralDof[] dofs, double[] values)
	{
		for (int i = 0; i < dofs.Length; i++)
		{
			constraints.Add(new NodalDisplacement(node, dofs[i], values[i]));
		}
	}
    
	
	public static void AssignConvectionDiffusionDirichletBCsToModel(Model model,
	List<(BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> bcs,
	double tolerance)
	{
		var constraints = new List<INodalConvectionDiffusionDirichletBoundaryCondition>();
		var emptyloads = new List<INodalConvectionDiffusionNeumannBoundaryCondition>();
		
		foreach (var bcSet in bcs)
		{
			var bcFace = bcSet.Item1;
			var bcDofs = bcSet.Item2;
			var faceCoords = bcSet.Item3;
			var bcValues = bcSet.Item4;
			
			foreach (var node in model.NodesDictionary)
			{
				var nodalCoords = new double[] { node.Value.X, node.Value.Y, node.Value.Z };
				
				if (bcFace == BoundaryConditionCase.RightDirichlet || bcFace == BoundaryConditionCase.LeftDirichlet
				    && Math.Abs(nodalCoords[0] - faceCoords[0][0]) < tolerance)
				{
					ApplyConvectionDiffusionDirichletBCToNode(constraints, node.Value, bcDofs, bcValues);
				}
				
				if (bcFace == BoundaryConditionCase.FrontDirichlet || bcFace == BoundaryConditionCase.BackDirichlet &&
				    Math.Abs(nodalCoords[1] - faceCoords[0][1]) < tolerance)
				{
					ApplyConvectionDiffusionDirichletBCToNode(constraints, node.Value, bcDofs, bcValues);
				}
				
				if (bcFace == BoundaryConditionCase.TopDirichlet || bcFace == BoundaryConditionCase.BottomDirichlet &&
				    Math.Abs(nodalCoords[2] - faceCoords[0][2]) < tolerance)
				{
					ApplyConvectionDiffusionDirichletBCToNode(constraints, node.Value, bcDofs, bcValues);
				}
			}
		}
		model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(constraints, emptyloads));
	}
	
	private static void ApplyConvectionDiffusionDirichletBCToNode(List<INodalConvectionDiffusionDirichletBoundaryCondition> constraints,
		INode node, ConvectionDiffusionDof[] dofs, double[] values)
	{
		constraints.Add(new NodalUnknownVariable(node, dofs[0], values[0]));
	}
}
