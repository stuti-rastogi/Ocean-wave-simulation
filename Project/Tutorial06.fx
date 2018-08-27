//--------------------------------------------------------------------------------------
// File: Tutorial06.fx
//
// Copyright (c) Microsoft Corporation. All rights reserved.
//--------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------
// Constant Buffer Variables
//--------------------------------------------------------------------------------------
cbuffer ConstantBuffer : register( b0 )
{
	matrix World;
	matrix View;
	matrix Projection;
	float4 vLightDir[2];
	float4 vLightColor[2];
	float4 vOutputColor;
}

//Light Structure and Variables//

struct DirectionalLight
{
	float4 color;
	float3 dir;
};

struct Material
{
	float Ka, Kd, Ks, specPow;
};

//--------------------------------------------------------------------------------------
struct VS_INPUT
{
    float4 Pos : POSITION;
    float3 Norm : NORMAL;
};

struct PS_INPUT
{
    float4 Pos : SV_POSITION;
    float3 Norm : TEXCOORD0;
	float4 worldPosition: POSITION1;
};


//--------------------------------------------------------------------------------------
// Vertex Shader
//--------------------------------------------------------------------------------------
PS_INPUT VS( VS_INPUT input )
{
    PS_INPUT output = (PS_INPUT)0;
    output.Pos = mul( input.Pos, World );
    output.Pos = mul( output.Pos, View );
    output.Pos = mul( output.Pos, Projection );
    output.Norm = mul( input.Norm, (float3x3)World );

	output.worldPosition = mul(input.Pos, World);
    
    return output;
}


float4 calcPhong(float4 LColor, float3 N, float3 L, float3 V, float3 R)
{
	float Ka = 0.5f;
	float Kd = 0.6f;
	float Ks = 1.f;
	int specPower = 2;
	float4 ambientLight;
	ambientLight = float4(0.f, 0.3f, 0.9f, 1.0f);
	float4 Ia = Ka * ambientLight;
	float4 Id = Kd * saturate(dot(N, L));
	float4 Is = Ks * pow(saturate(dot(R, V)), specPower);

	//return float4(Id);
	return float4(Ia + ((Id + Is)) * LColor);
}

//--------------------------------------------------------------------------------------
// Pixel Shader
//--------------------------------------------------------------------------------------
float4 PS( PS_INPUT input) : SV_Target
{

	DirectionalLight light;
	light.dir = normalize(float3(1,-1,0));
	light.color = normalize(float4(1.f,1.f,1.f,1.0f));

    float4 finalColor = 0;
    
	input.Norm = normalize(input.Norm);
	float3 V = normalize(float3(-21.f, 10.f, -21.f) - (float3)input.worldPosition.xyz);

	float3 R = reflect(light.dir, input.Norm);

	float4 I = calcPhong(light.color, input.Norm, -light.dir, V, R);

    ////do NdotL lighting for 2 lights
    //for(int i=0; i<2; i++)
    //{
    //    finalColor += saturate( dot( (float3)vLightDir[i],input.Norm) * vLightColor[i] );
    //}
    //finalColor.a = 1;
    return I;
	
}


//--------------------------------------------------------------------------------------
// PSSolid - render a solid color
//--------------------------------------------------------------------------------------
float4 PSSolid( PS_INPUT input) : SV_Target
{
    return vOutputColor;
}
