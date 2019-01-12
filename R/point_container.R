pointContainer <- setRefClass("pointContainer",
                      fields = list(p = "numeric",
                                    L = "function",
                                    input.dimension = "numeric",
                                    epsilon.shift.input = "numeric",
                                    EdgeDistance = "function"))
