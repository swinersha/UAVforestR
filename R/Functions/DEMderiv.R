DEMderiv<-function(data,attr,method){
  #Title: Function for computing terrain attributes from a raster
  #Author: Fabio Veronesi
  #License: Creative Commons - Attribution-NonCommercial (CC BY-NC) - http://creativecommons.org/
  
  cellValue=res(data)[1]
  
  evans<-function(neighb){  #method proposed by Evans (1980) as it is described in Florisky (1998)
    #Node numbering for the 3x3 window
    #Z1    Z2    Z3
    #Z4    Z5    Z6
    #Z7    Z8    Z9
    z1<-neighb[1]
    z2<-neighb[2]
    z3<-neighb[3]
    z4<-neighb[4]
    z6<-neighb[6]
    z7<-neighb[7]
    z8<-neighb[8]
    z9<-neighb[9]
    z5<-neighb[5]
    r=(z1+z3+z4+z6+z7+z9-(2*(z2+z5+z8)))/(3*(cellValue^2))
    t=(z1+z2+z3+z7+z8+z9-(2*(z4+z5+z6)))/(3*(cellValue^2))
    s=(z3+z7-z1-z9)/(4*(cellValue^2))
    p=(z3+z6+z9-z1-z4-z7)/(6*cellValue)
    q=(z1+z2+z3-z7-z8-z9)/(6*cellValue)
    if(paste(attr)=="slope"){result= atan(sqrt(p^2+q^2)) }                 
    else{
      if(paste(attr)=="aspect"){result= 180-atan2(q,p)+90*(p/abs(p)) }                         
      else{
        if(paste(attr)=="plan.curvature"){result= -(q^2*r-2*p*q*s+p^2*t)/((p^2+q^2)*sqrt(1+p^2+q^2))}      
        else{
          if(paste(attr)=="prof.curvature"){result= -(p^2*r+2*p*q*s+q^2*t)/((p^2+q^2)*sqrt(1+p^2+q^2)^3)}  
        }}}
    c(result)
  }
  
  shary<-function(neighb){  #method proposed by Shary (1995) as it is described in Florisky (1998)
    #Node numbering for the 3x3 window
    #Z1    Z2    Z3
    #Z4    Z5    Z6
    #Z7    Z8    Z9
    z1<-neighb[1]
    z2<-neighb[2]
    z3<-neighb[3]
    z4<-neighb[4]
    z6<-neighb[6]
    z7<-neighb[7]
    z8<-neighb[8]
    z9<-neighb[9]
    z5<-neighb[5]
    r=(z1+z3+z7+z9+3*(z4+z6)-2*(z2+3*z5+z8))/(5*(cellValue^2))
    t=(z1+z3+z7+z9+3*(z2+z8)-2*(z4+3*z5+z6))/(5*(cellValue^2))
    s=(z3+z7-z1-z9)/(4*(cellValue^2))
    p=(z3+z6+z9-z1-z4-z7)/(6*cellValue)
    q=(z1+z2+z3-z7-z8-z9)/(6*cellValue)             
    if(paste(attr)=="slope"){result= atan(sqrt(p^2+q^2)) }                 
    else{
      if(paste(attr)=="aspect"){result= 180-atan2(q,p)+90*(p/abs(p)) }                         
      else{
        if(paste(attr)=="plan.curvature"){result= -(q^2*r-2*p*q*s+p^2*t)/((p^2+q^2)*sqrt(1+p^2+q^2)) }       
        else{
          if(paste(attr)=="prof.curvature"){result= -(p^2*r+2*p*q*s+q^2*t)/((p^2+q^2)*sqrt(1+p^2+q^2)^3) }  
        }}}
    c(result)
  }
  
  
  zev.tho<-function(neighb){   #method proposed by Zevenbergen and Thorne (1987) as it is described in Florisky (1998)
    #Node numbering for the 3x3 window
    #Z1    Z2    Z3
    #Z4    Z5    Z6
    #Z7    Z8    Z9
    z1<-neighb[1]
    z2<-neighb[2]
    z3<-neighb[3]
    z4<-neighb[4]
    z6<-neighb[6]
    z7<-neighb[7]
    z8<-neighb[8]
    z9<-neighb[9]
    z5<-neighb[5]
    p=(z6-z4)/(2*cellValue)
    q=(z2-z8)/(2*cellValue)
    r=(z4+z6-2*z5)/(2*(cellValue^2))
    s=(z3+z7-z1-z9)/(4*(cellValue^2))
    t=(z2+z8-2*z5)/(2*(cellValue^2))
    if(paste(attr)=="slope"){result= atan(sqrt(p^2+q^2)) }                 
    else{
      if(paste(attr)=="aspect"){result= 180-atan2(q,p)+90*(p/abs(p)) }                         
      else{
        if(paste(attr)=="plan.curvature"){result= -(q^2*r-2*p*q*s+p^2*t)/((p^2+q^2)*sqrt(1+p^2+q^2))  }       
        else{
          if(paste(attr)=="prof.curvature"){result= -(p^2*r+2*p*q*s+q^2*t)/((p^2+q^2)*sqrt(1+p^2+q^2)^3) }  
        }}}
    c(result)
  }
  
  
  
  moore<-function(neighb){   #method proposed by Moore et al. (1993) as it is described in Florisky (1998)
    #Node numbering for the 3x3 window
    #Z1    Z2    Z3
    #Z4    Z5    Z6
    #Z7    Z8    Z9
    z1<-neighb[1]
    z2<-neighb[2]
    z3<-neighb[3]
    z4<-neighb[4]
    z6<-neighb[6]
    z7<-neighb[7]
    z8<-neighb[8]
    z9<-neighb[9]
    z5<-neighb[5]
    p=(z6-z4)/(2*cellValue)
    q=(z2-z8)/(2*cellValue)
    r=(z4+z6-2*z5)/(cellValue^2)
    s=(z3+z7-z1-z9)/(4*(cellValue^2))
    t=(z2+z8-2*z5)/(cellValue^2)
    if(paste(attr)=="slope"){result= atan(sqrt(p^2+q^2)) }                 
    else{
      if(paste(attr)=="aspect"){result= 180-atan2(q,p)+90*(p/abs(p)) }                         
      else{
        if(paste(attr)=="plan.curvature"){result= -(q^2*r-2*p*q*s+p^2*t)/((p^2+q^2)*sqrt(1+p^2+q^2)) }       
        else{
          if(paste(attr)=="prof.curvature"){result= -(p^2*r+2*p*q*s+q^2*t)/((p^2+q^2)*sqrt(1+p^2+q^2)^3) }  
        }}}
    c(result)
  }
  
  
  
  if(paste(method)=="evans"){crop(focal(data,w=matrix(1,3,3),fun=evans,pad=T,padValue=0),y=c(data@extent@xmin+cellValue,data@extent@xmax-cellValue,data@extent@ymin+cellValue,data@extent@ymax-cellValue))}
  else{
    if(paste(method)=="zev.tho"){crop(focal(data,w=matrix(1,3,3),fun=zev.tho,pad=T,padValue=0),y=c(data@extent@xmin+cellValue,data@extent@xmax-cellValue,data@extent@ymin+cellValue,data@extent@ymax-cellValue))}
    else{
      if(paste(method)=="shary"){crop(focal(data,w=matrix(1,3,3),fun=shary,pad=T,padValue=0),y=c(data@extent@xmin+cellValue,data@extent@xmax-cellValue,data@extent@ymin+cellValue,data@extent@ymax-cellValue))}
      else{
        if(paste(method)=="moore"){crop(focal(data,w=matrix(1,3,3),fun=moore,pad=T,padValue=0),y=c(data@extent@xmin+cellValue,data@extent@xmax-cellValue,data@extent@ymin+cellValue,data@extent@ymax-cellValue))}
      }}}
  
  
  #REFERENCES
  #Evans, I. S. (1980). An Integrated System of Terrain Analysis and Slope Mapping. Zeitschrift für Geomorphologie, Suppl. Bd. 36, 274-295.
  #Florinsky, I. V. (1998). Accuracy of Local Topographic Variables Derived from Digital Elevation Models. International Journal of Geographical Informaation Science, 12:1, 47-62.
  #Shary, P. A. (1991). The Second Derivative Topographic Method. In: The Geometry of the Earth Surface Structures, edited by Stepanov, I. N., 15-29 (in Russian).
  #Wilson, J. P. & Gallant, J. C. (2000). Terrain Analysis - Principles and Applications. Wiley.
  #Zevenbergen, L. W. & Thorne, C. R. (1987). Quantitative Analysis of Land Surface Topography. Earth Surface Processes and Landforms, vol. 12, 47-56.
  #Moore, I. D.; Gessler, P. E.; Nielsen, G. A. & Paterson, G. A. (1993). Soil Attribute Prediction Using Terrain Analysis. Soil Science Society of America Journal, vol. 57, 443-452.
}