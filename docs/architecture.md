# Role 5 아키텍처

```
[Component 잔차 + 포트 연결식] -->
[system/assembly.py] -- Global F(t,y,ydot)=0 -->
[solvers/dassl.py] -- 가변 간격 BDF 1–5차 + Newton --> y(t), ydot(t)
```

- `interfaces/`: 모든 개발자의 변수 등록 및 단위/Port 규약
- `system/`: Component 잔차 뒤에 포트 연결식을 조립하고 전역 미지수와 잔차의 개수를 확인
- `solvers/`: index-1 DAE 일치 초기화, 정상상태 계산, 시간적분
- `examples/`: 실제 물리 모델에 앞서 구조를 테스트하는 별도의 가상 DAE
- `models/`, `physics/`: 추후 담당 팀원이 구현할 부분 (현재는 명시적 placeholder)

## 시간적분

첫 스텝은 후진 오일러와 절반 스텝을 사용해 국소오차를 추정한다. 이후에는 수용된
시각과 상태로 보간 다항식을 만들고, 차수 `k`의 BDF에서 Newton 보정에 곱하는
계수를 `(1 + 1/2 + … + 1/k) / h`로 둔다. 가중 RMS 오차가 1 이하인 스텝을
수용한 뒤 인접 차수의 오차 추정치를 비교해 다음 간격과 차수를 정한다. 수용된
이력만 다음 스텝에 사용한다. [Petzold의 DASSL 설명](https://cse.cs.ucsb.edu/sites/default/files/publications/DASSL.pdf)

Python API는 NumPy 밀집 선형계와 수용 스텝 출력을 제공한다. 수치 차분 Jacobian이
기본값이며, 사용자가 결합 Jacobian을 전달할 수 있다.
