using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.DrugDeliveryModel.Tests.Commons;

public class BoundaryConditionsUtility
{
	public enum BoundaryConditionCase
	{
		LeftDirichlet,
		LeftPointForce,
		LeftDistributedForce,
		RightDirichlet,
		RightPointForce,
		RightDistributedForce,
		FrontDirichlet,
		FrontPointForce,
		FrontDistributedForce,
		BackDirichlet,
		BackPointForce,
		BackDistributedForce,
		BottomDirichlet,
		BottomPointForce,
		BottomDistributedForce,
		TopDirichlet,
		TopPointForce,
		TopDistributedForce
	}

	public static void AssignStructuralDirichletBCs(Model model,
		List<(BoundaryConditionCase, StructuralDof[], double[][], double[])> bcs,
		double tolerance)
	{

		var constraints = new List<INodalDisplacementBoundaryCondition>();
		var emptyloads = new List<INodalLoadBoundaryCondition>();


		foreach (var node in model.NodesDictionary)
		{
			var nodalCoords = new double[] { node.Value.X, node.Value.Y, node.Value.Z };

			foreach (var bcSet in bcs)
			{
				var bcFace = bcSet.Item1;
				var bcDofs = bcSet.Item2;
				var faceCoords = bcSet.Item3;
				var bcValues = bcSet.Item4;

				if (bcFace == BoundaryConditionCase.RightDirichlet || bcFace == BoundaryConditionCase.LeftDirichlet)
				{
					if (nodalCoords[0] > faceCoords[0][0] - tolerance)
					{
						ApplyStructuralDirichletBCToNode(constraints, node.Value, bcDofs, bcValues);
					}
				}
				else if (bcFace == BoundaryConditionCase.FrontDirichlet ||
				         bcFace == BoundaryConditionCase.BackDirichlet)
				{
					if (nodalCoords[1] > faceCoords[0][1] - tolerance)
					{
						ApplyStructuralDirichletBCToNode(constraints, node.Value, bcDofs, bcValues);
					}
				}
				else if (bcFace == BoundaryConditionCase.TopDirichlet ||
				         bcFace == BoundaryConditionCase.BottomDirichlet)
				{
					if (nodalCoords[2] > faceCoords[0][2] - tolerance)
					{
						ApplyStructuralDirichletBCToNode(constraints, node.Value, bcDofs, bcValues);
					}
				}
				else
				{
					throw new Exception(
						"BoundaryConditionCase should be LeftDirichlet, RightDirichlet, FrontDirichlet, BackDirichlet, TopDirichlet or BottomDirichlet");
				}
			}

			model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, emptyloads));
		}
	}

	public static void AssignStructuralNeumannBCs(Model model, List<(BoundaryConditionCase, StructuralDof[], double[][], double[])> bcs,
	    double tolerance)
    {
	    var emptyConstraints = new List<INodalDisplacementBoundaryCondition>();
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

                if (loadCase == BoundaryConditionCase.RightPointForce ||
                    loadCase == BoundaryConditionCase.LeftPointForce)
                {
                    for (int i = 0; i <= loadCoords.Length - 1; i++)
                    {
	                    if (nodalCoords[0] > loadCoords[i][0] - tolerance)
	                    {
		                    ApplyStructuralNeumannBCToNode(loads, node.Value, bcDofs, bcValues);
	                    }
                    }
                }

                else if (loadCase == BoundaryConditionCase.FrontPointForce ||
                         loadCase == BoundaryConditionCase.BackPointForce)
                {
                    for (int i = 0; i <= loadCoords.Length - 1; i++)
                    {
	                    if (nodalCoords[1] > loadCoords[i][1] - tolerance)
	                    {
		                    ApplyStructuralNeumannBCToNode(loads, node.Value, bcDofs, bcValues);
	                    }
                    }
                }
                
                else if (loadCase == BoundaryConditionCase.TopPointForce ||
                         loadCase == BoundaryConditionCase.BottomPointForce)
                {
                    for (int i = 0; i <= loadCoords.Length - 1; i++)
                    {
	                    if (nodalCoords[2] > loadCoords[i][2] - tolerance)
	                    {
		                    ApplyStructuralNeumannBCToNode(loads, node.Value, bcDofs, bcValues);
	                    }
                    }
                }
                
                else if (loadCase == BoundaryConditionCase.RightDistributedForce ||
                         loadCase == BoundaryConditionCase.LeftDistributedForce)
                {
                    if (nodalCoords[0] > loadCoords[0][0] - tolerance)
                    {
	                    ApplyStructuralNeumannBCToNode(loads, node.Value, bcDofs, bcValues);
                    }
                }
                
                else if (loadCase == BoundaryConditionCase.FrontDistributedForce ||
                         loadCase == BoundaryConditionCase.BackDistributedForce)
                {
                    if (nodalCoords[1] > loadCoords[0][1] - tolerance)
                    {
	                    ApplyStructuralNeumannBCToNode(loads, node.Value, bcDofs, bcValues);
                    }
                }
                
                else if (loadCase == BoundaryConditionCase.TopDistributedForce ||
                         loadCase == BoundaryConditionCase.BottomDistributedForce)
                {
                    if (nodalCoords[2] > loadCoords[0][2] - tolerance)
                    {
	                    ApplyStructuralNeumannBCToNode(loads, node.Value, bcDofs, bcValues);
                    }
                }
                
                else
                {
                    throw new Exception 
	                    ("For Neumann Boundary Conditions, the face must be: Left, Right, Front, Back, Bottom, Top");
                }

                ;
            }
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(emptyConstraints, loads));
        }
    }
    
    public static void ApplyStructuralNeumannBCToNode(List<INodalLoadBoundaryCondition> constraints,
	    INode node, StructuralDof[] dofs, double[] values)
    {
	    for (int i = 0; i < dofs.Length; i++)
	    {
		    constraints.Add(new NodalLoad(node, dofs[i], values[i]));
	    }
    }

	public static void ApplyStructuralDirichletBCToNode(List<INodalDisplacementBoundaryCondition> constraints,
		INode node, StructuralDof[] dofs, double[] values)
	{
		for (int i = 0; i < dofs.Length; i++)
		{
			constraints.Add(new NodalDisplacement(node, dofs[i], values[i]));
		}
	}
}
