
answer1arg() {
  answer=$(../fino $1 $2)
  
  if [ "${answer}" = "$3" ]; then
    echo "ok"
    level=0
  else
    echo "wrong"
    level=1
  fi

  return ${level}
}


answer2argssorthead1() {
  answer=$(../fino $1 $2 $3 | sort -rg | head -n1)
  
  if [ "${answer}" = "$4" ]; then
    echo "ok"
    level=0
  else
    echo "wrong"
    level=1
  fi

  return ${level}
}
